/*
            HYDRA
        Version 2.0
    last edited: 25.01.2021
    Michaela Schneeberger,
    Dr. Philipp Bahlke,
    Prof. Dr. Carmen Herrmann
    Universitaet Hamburg
*/

#define _USE_MATH_DEFINES

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <sstream>
#include <complex>
#include <mpi.h>
#include <omp.h>
#include "InputOutput.hpp"

#define USAGE\
    std::cerr << "Usage: " << argv[0] << "Execute without arguments in directory with input \
        (eigenvalues, coef, hamiltonian, overlap, basis_array, delta.global.in)"; \
    exit(EXIT_FAILURE);

/**
 * Build and generate energy grid
 */
Eigen::ArrayXcd build_Egrid(int points, int beta, float step, float smear, int grid_type)
{
    Eigen::ArrayXcd E = Eigen::ArrayXcd::Zero(points);

    if (grid_type == 2) {
        for (int i = 0;  i < points; i++) {
            std::complex<double> value((double) 0.0, (double) (((2*(i-1))+1)*M_PI)/beta);
            E[i] = value;
        }
    }
    else if (grid_type == 1) {
        for (int i = 0; i < points; i++ ) {
            std::complex<double> value((double) (+1*(i+1)*step)-((points/2)*step), (double) smear);
            E[i] = value;
        }
    }
    return E;
}

//From local Greens function gf compute delta_ii = gf_ii - w, 
Eigen::ArrayXcd calc_delta(Eigen::MatrixXcd gf, int ncorrorbs, int xbf1, Eigen::MatrixXd over, double homo, \
                            std::complex<double> E, int leh)
{
    Eigen::MatrixXd S_core;
    Eigen::MatrixXcd gf_sub, gf_sub_inv, Iden, delta_mat;
    Eigen::ArrayXcd delta = Eigen::ArrayXcd::Zero(ncorrorbs);

    delta_mat = Eigen::MatrixXcd::Zero(ncorrorbs, ncorrorbs);
    S_core = over.block(xbf1 - 1, xbf1 - 1, ncorrorbs, ncorrorbs);

    if (leh == 1)
        gf_sub = gf;
    else
        gf_sub = gf.block(xbf1 - 1, xbf1 - 1, ncorrorbs, ncorrorbs);

    //test if inversion was succesfull
    gf_sub_inv = gf_sub.inverse();
    Iden = gf_sub * gf_sub_inv;

    for (int i = 0; i < ncorrorbs; i++) {
        if (1.000 > Iden(i, i).real() && 0.999 < Iden(i, i).real() && 0.00001 > std::abs(Iden(i, i).imag())) {
            std::complex<double> c(1.0, 0.0);
            Iden(i, i) = c;
        }
    }

    for (int i = 0; i < ncorrorbs; i++) {
        if (Iden(i, i).real() != 1.0 && std::abs(Iden(i, i).imag()) > 0.0001) {
            std::cerr << "Inversion error at energy " << E << " for lehmann " << leh << " : (GF*GF^-1)_" \
            << i << i << " = " << Iden(i, i);
            exit(EXIT_FAILURE);
        }
    }

    //calulate delta
    if (leh == 1) {
        for (int i = 0; i < ncorrorbs; i++) {
            for (int j = 0; j < ncorrorbs; j++) {
                delta_mat(i, j) = -(gf_sub_inv(i, j) - (std::complex < double > )(E - homo));
            }
        }
    } else {
        for (int i = 0; i < ncorrorbs; i++) {
            for (int j = 0; j < ncorrorbs; j++) {
                delta_mat(i, j) = -(gf_sub_inv(i, j) - (std::complex < double > )((E - homo) * S_core(i, j)));
            }
        }
    }

    for (int i = 0; i < ncorrorbs; i++) {
        delta(i) = delta_mat(i, i);
    }

    return delta;
}

//get local Greens function gf by inversion of full Greens function
Eigen::MatrixXcd make_direct_GF(Eigen::MatrixXd &Ham_mod, Eigen::MatrixXd &Over_mod, int nbas, double homo, \
                                std::complex<double> E)
{
    Eigen::MatrixXcd gf = Eigen::MatrixXcd::Zero(nbas, nbas);

    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            gf(i, j) = ((E - homo) * Over_mod(i, j) - Ham_mod(i, j));
        }
    }

    gf = gf.inverse();
    return gf;

}

//Loewdin orthogonalization
Eigen::MatrixXd LoewdinOrtho(Eigen::MatrixXd &M, int n)
{
    Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> es(M);
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; i++) {
        double x = es.eigenvalues()[i];
        double y = std::sqrt(x);
        X(i, i) = 1 / y;
    }

    X = es.eigenvectors() * X * es.eigenvectors().transpose();
    return X;
}

//print orbitals in molden compatible format
void print_orbitals(Eigen::MatrixXd mat, Eigen::VectorXcd eig, int norbs, double homo, int output, int xbf1, int nbas)
{
    std::ofstream molden;
    if (output == 1) {
        molden.open("molden.impurity");
    } else {
        molden.open("molden.bath");
    }
    for (int i = 0; i < norbs; ++i) {
        molden << "Sym= " << i + 1 << " a" << std::endl;
        molden << "Ene= " << std::real(eig(i) - homo) << std::endl;
        molden << "Spin= Alpha" << std::endl;
        molden << "Occup= 0.000000" << std::endl;
        if (output == 1) {
            for (int j = 0; j < xbf1; ++j) {
                molden << "        " << j + 1 << "    " << "0.0" << std::endl;
            }
            for (int j = xbf1; j < (xbf1+norbs); ++j) {
                molden << "        " << j + 1 << "    " << mat(j-xbf1, i) << std::endl;
            }
            for (int j = (xbf1+norbs); j < nbas; ++j) {
                molden << "        " << j + 1 << "    " << "0.0" << std::endl;
            }
        } else {
            for (int j = 0; j < xbf1; ++j) {
                molden << "        " << j + 1 << "    " << mat(j, i) << std::endl;
            }
            for (int j = xbf1; j < (xbf1-norbs+nbas); ++j) {
                molden << "        " << j + 1 << "    " << "0.0" << std::endl;
            }
            for (int j = (xbf1-norbs+nbas); j < nbas; ++j) {
                molden << "        " << j + 1 << "    " << mat(j+norbs-nbas, i) << std::endl;
            }
        }
    }
    molden.close();
}

//diagonalize impurity orbitals, return diagonalization matrix
Eigen::MatrixXd diag_corr_orbs(Eigen::MatrixXd &ham, int xbf1, int norb, double homo, int nbas, int mpi) {

    Eigen::MatrixXd trans(nbas, nbas);
    trans.setIdentity();

    Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> es(ham.block(xbf1 - 1, xbf1 - 1, norb, norb));
    trans.block(xbf1 - 1, xbf1 - 1, norb, norb) = es.eigenvectors().real();

    if (mpi == 0) {
        std::cout << "Rotated basis" << std::endl;
        if (norb == 5) {
            std::cout << "z^2		xz		yz		xy		x^2-y^2" << std::endl;
        } else if (norb == 3) {
            std::cout << "px        py      pz" << std::endl;
        }
        std::cout << trans.block(xbf1 - 1, xbf1 - 1, norb, norb);
    }

    if (mpi == 0) {
        print_orbitals(trans.block(xbf1 - 1, xbf1 - 1, norb, norb), es.eigenvalues(), norb, homo, 1, xbf1, nbas);
    }

    return trans;
}

//diagonalize bath orbitals, return diagonalization matrix
Eigen::MatrixXd diag_bath_orbs(int xbf1, int xbf2, int nbas, int ncorrorbs, Eigen::MatrixXd &Ham_ortho, double homo)
{
    Eigen::MatrixXd coup_Trans, block_mat;
    coup_Trans = Eigen::MatrixXd::Identity(nbas, nbas);

    block_mat = Eigen::MatrixXd::Zero(nbas - ncorrorbs, nbas - ncorrorbs);
    block_mat.block(0, 0, xbf1 - 1, xbf1 - 1) = Ham_ortho.block(0, 0, xbf1 - 1, xbf1 - 1);
    block_mat.block(0, xbf1 - 1, xbf1 - 1, nbas - xbf2) = Ham_ortho.block(0, xbf2, xbf1 - 1, nbas - xbf2);
    block_mat.block(xbf1 - 1, 0, nbas - xbf2, xbf1 - 1) = Ham_ortho.block(xbf2, 0, nbas - xbf2, xbf1 - 1);
    block_mat.block(xbf1 - 1, xbf1 - 1, nbas - xbf2, nbas - xbf2) = Ham_ortho.block(xbf2, xbf2, nbas - xbf2,
                                                                                    nbas - xbf2);
    Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> es(block_mat);
    block_mat = es.eigenvectors().real();

    print_orbitals(block_mat, es.eigenvalues(), nbas-ncorrorbs, homo, 2, xbf1, nbas);

    coup_Trans.block(0, 0, xbf1 - 1, xbf1 - 1) = block_mat.block(0, 0, xbf1 - 1, xbf1 - 1);
    coup_Trans.block(0, xbf2, xbf1 - 1, nbas - xbf2) = block_mat.block(0, xbf1 - 1, xbf1 - 1, nbas - xbf2);
    coup_Trans.block(xbf2, 0, nbas - xbf2, xbf1 - 1) = block_mat.block(xbf1 - 1, 0, nbas - xbf2, xbf1 - 1);
    coup_Trans.block(xbf2, xbf2, nbas - xbf2, nbas - xbf2) = block_mat.block(xbf1 - 1, xbf1 - 1, nbas - xbf2,
                                                                             nbas - xbf2);
    return coup_Trans;
}

//build projector
Eigen::MatrixXd buildproj(Eigen::MatrixXd &Over, Eigen::MatrixXd Coef, int nbas, Eigen::MatrixXd &Eig, \
                            int diag_param, int xbf1, int ncorrorbs, Eigen::MatrixXd &X, \
                            int mpi_rank) {
    Eigen::MatrixXd P(nbas, nbas), temp(nbas, nbas), P_C(nbas, nbas);

    P_C.setIdentity();

    P = X * Over * Coef.transpose();

    temp = P.transpose() * Eig * P;

    if (diag_param == 1) {
        temp = diag_corr_orbs(temp, xbf1, ncorrorbs, 0.0, nbas, mpi_rank);
        P = temp.transpose() * P;
    }

    P_C.block(xbf1 - 1, 0, ncorrorbs, nbas) = P.block(xbf1 - 1, 0, ncorrorbs, nbas);

    return P_C;
}

//build full Greens function GF in Lehmann representation 
Eigen::MatrixXcd make_lehmann(Eigen::MatrixXd Eig, int nbas, std::complex<double> E, \
                                double homo) {
    Eigen::MatrixXcd GF(nbas, nbas);
    for (int i = 0; i < nbas; ++i) {
        GF(i, i) = 1.0 / ((E - homo) - Eig(i, i));
    }
    return GF;
}

//get local Greens functin gf by projection on full Greens function GF in Lehmann representation
Eigen::MatrixXcd projection(Eigen::MatrixXd P_C, Eigen::MatrixXcd &GF, int nbas, \
                            int xbf1, int ncorrorbs) {
    Eigen::MatrixXcd gf(nbas, nbas);
    gf = P_C.transpose() * GF * P_C;

    return gf.block(xbf1 - 1, xbf1 - 1, ncorrorbs, ncorrorbs);
}


int main(int argc, char *argv[]) {

    if (argc != 1) {
        USAGE;
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    //set up mpi and omp
    MPI::Init(argc, argv);
    MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
    int mpi_rank = 0;
    int mpi_size = 0;
    try {
        mpi_rank = MPI::COMM_WORLD.Get_rank();
        mpi_size = MPI::COMM_WORLD.Get_size();
    }
    catch (MPI::Exception e) {
        std::cout << "MPI_Error: " << e.Get_error_code() << " - " << e.Get_error_string() << std::endl;
    }
    int omp_size = omp_get_max_threads();

    //Read user input from file delta.global.in
    std::fstream infile;
    infile.open("delta.global.in");
    std::vector<float> Param = read_input(infile, mpi_rank);
    infile.close();

    int nbas = Param[0]; //number of basis functions in whole system
    int xbf1 = Param[1]; //position of first basis function of corr. system 
    int xbf2 = Param[2]; //position of last basis function of corr. system
    int grid_param = Param[3]; //energy grid type: 0:real, 1:matsubara
    int beta = Param[4]; //beta for matsubara energy
    int points = Param[5]; //number of energy points to be used
    float step = Param[6]; //step size of energy grid
    float smear = Param[7]; //smear of real energy grid
    double homo = Param[8]; //energy of HOMO orbital, will be shiftet to w=0
    int unit = Param[9]; //0: eV, 1: Hartree
    int diag_param = Param[10]; //diagonalize corr. system? 1:yes 0:no
    int gf_param = Param[11]; //which method? 0:bad 1:Green's 2:projector(not fully implemented)
    int coup_switch = Param[12]; //calculate hybridization of subsystem specified in file basis_array, 1:yes 0:no

    //even out number of energy points between mpi processes
    int mpipoints = points / mpi_size;
    points = (mpipoints + 1) * mpi_size;
    int my_points = points / mpi_size;

    int ncorrorbs = xbf2 - xbf1 + 1; //number of basis functions in corr. subsystem

    //Print information of current run, see InputOutput.cpp
    if (mpi_rank == 0) {
        print_output(nbas, xbf1, xbf2, grid_param, beta, points, step, smear, homo, unit, diag_param, gf_param, \
         coup_switch, my_points, mpi_size, omp_size);
    }

    //Read subsystem basis functions from file basis_array to calculate decomposition of hybrid function
    Eigen::ArrayXd b_arr(nbas);
    if (coup_switch == 1) {
        infile.open("basis_array");
        if (mpi_rank == 0) {
            std::cout << "Reading basis_array\n";
        }
        b_arr = read_intvector(infile, nbas, mpi_rank);
        infile.close();
    }

    //Read coefficient matrix from file coef in variable Coef
    Eigen::MatrixXd Coef;
    infile.open("coef");
    if (mpi_rank == 0) {
        std::cout << "Reading coef matrix" << std::endl;
    }
    Coef = read_denmat(infile, nbas, mpi_rank, "coef");

    //Read hamiltonian matrix from file hamiltonian in variable Ham
    Eigen::MatrixXd Ham;
    infile.open("hamiltonian");
    if (mpi_rank == 0) {
        std::cout << "Reading Hamiltonian matrix" << std::endl;
    }
    Ham = read_denmat(infile, nbas, mpi_rank, "hamiltonian");

    //Read overlap matrix from file overlap in variable Over
    Eigen::MatrixXd Over;
    infile.open("overlap");
    if (mpi_rank == 0) {
        std::cout << "Reading overlap matrix" << std::endl;
    }
    Over = read_denmat(infile, nbas, mpi_rank, "overlap");

    //Read eigenvalues matrix from file eigenvalues in variable Eig
    Eigen::MatrixXd Eig;
    infile.open("eigenvalues");
    if (mpi_rank == 0) {
        std::cout << "Reading eigenvalues" << std::endl;
    }
    Eig = read_denmat(infile, nbas, mpi_rank, "eigenvalues");

    //store energy grid in variable Earray
    Eigen::ArrayXcd Earray = build_Egrid(points, beta, step, smear, grid_param);

    //convert Energy of hamiltonian and eigenvalues matrices if needed
    double e_conv;
    if (unit == 1)
        e_conv = 27.21138;

    Ham = Ham * e_conv;
    Eig = Eig * e_conv;

    if (mpi_rank == 0) {
        auto t4 = std::chrono::high_resolution_clock::now();
        auto duration3 = std::chrono::duration_cast<std::chrono::seconds>(t4 - t1).count();
        std::cout << "Time until matrices are read: " << duration3 << std::endl;
    }

    Eigen::MatrixXd Ham_ortho, X;
    X = Eigen::MatrixXd::Zero(nbas, nbas);

    //Build matrix X to Lowedin orthogonalize atom-centered basis functions
    Eigen::MatrixXd M = Over;
    X = LoewdinOrtho(M, nbas);
    Ham_ortho = X * Ham * X;

    //Build matrix T to diagonalize (and ortho.) correlated sub space if requested
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(nbas, nbas);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nbas, nbas);
    if (diag_param == 1) {
        D = diag_corr_orbs(Ham_ortho, xbf1, ncorrorbs, 0.0, nbas, mpi_rank);
        T = X * D;
    } else {
        T = X;
    }
    Eigen::MatrixXd Ham_ortho_dia = T.transpose() * Ham * T;
    Eigen::MatrixXd Over_ortho_dia = T.transpose() * Over * T;

    //calculate projector for Lehmann representation (does not work yet!)
    Eigen::MatrixXd P_C(nbas, nbas);
    if (gf_param == 2) {
        P_C = buildproj(Over, Coef, nbas, Eig, diag_param, xbf1, ncorrorbs, X, mpi_rank);
    }

    if (mpi_rank == 0) {
        auto t5 = std::chrono::high_resolution_clock::now();
        auto duration4 = std::chrono::duration_cast<std::chrono::seconds>(t5 - t1).count();
        std::cout << "Time used before delta: " << duration4 << std::endl;
    }

    if (mpi_rank == 0) {
        std::cout << "Calculate Greensfunctions and Delta\n";
    }

    //Projector or Greens function approach
    if (gf_param == 1 || gf_param == 2) {
        //Variables for each mpi process
        Eigen::MatrixXcd my_gf_proj_dia = Eigen::MatrixXcd::Zero(ncorrorbs, my_points);
        Eigen::MatrixXcd my_D_lehmann = Eigen::MatrixXcd::Zero(ncorrorbs, my_points);
        Eigen::MatrixXcd tmp_direct_GF = Eigen::MatrixXcd::Zero(nbas, nbas);
        Eigen::MatrixXcd tmp_proj_GF = Eigen::MatrixXcd::Zero(ncorrorbs, ncorrorbs);
        Eigen::ArrayXcd temp_direct_Lehmann = Eigen::ArrayXcd::Zero(ncorrorbs);

        Eigen::ArrayXd E_imp_orbs = Eigen::ArrayXd::Zero(ncorrorbs);
        for (int i = 0; i < ncorrorbs; i++)
            E_imp_orbs(i) = Ham_ortho_dia(xbf1 - 1 + i, xbf1 - 1 + i);

        //Variables for complete matrices (for mpi rank 0)
        Eigen::MatrixXcd D_lehmann;
        Eigen::MatrixXcd gf_proj_dia;
        if (mpi_rank == 0) {
            D_lehmann = Eigen::MatrixXcd::Zero(ncorrorbs, points);
            gf_proj_dia = Eigen::MatrixXcd::Zero(ncorrorbs, points);
        }


        //divide energy grid between mpi processes
        try {
            MPI_Scatter(D_lehmann.data(), ncorrorbs * my_points, MPI::DOUBLE_COMPLEX, my_D_lehmann.data(), \
            ncorrorbs * my_points, MPI::DOUBLE_COMPLEX, 0, MPI::COMM_WORLD);
        }
        catch (MPI::Exception fail) {
            std::cerr << "MPI_Error: " << fail.Get_error_code() << " - " << fail.Get_error_string() << std::endl;
        }

        try {
            MPI_Scatter(gf_proj_dia.data(), ncorrorbs * my_points, MPI::DOUBLE_COMPLEX, my_gf_proj_dia.data(), \
           ncorrorbs * my_points, MPI::DOUBLE_COMPLEX, 0, MPI::COMM_WORLD);
        }
        catch (MPI::Exception fail) {
            std::cerr << "MPI_Error: " << fail.Get_error_code() << " - " << fail.Get_error_string() << std::endl;
        }

        //divide parts of energy grid between omp threads 
#pragma omp parallel
        {
#pragma omp for private(tmp_direct_GF, tmp_proj_GF, temp_direct_Lehmann) schedule(static, 1)
            for (int k = mpi_rank * my_points; k < ((mpi_rank + 1) * my_points); k++) {

                if (gf_param == 2) {
                    tmp_direct_GF = make_lehmann(Eig, nbas, Earray(k), -homo);
                    tmp_proj_GF = projection(P_C, tmp_direct_GF, nbas, xbf1, ncorrorbs);

                } else if (gf_param == 1) {
                    tmp_direct_GF = make_direct_GF(Ham_ortho_dia, Over_ortho_dia, nbas, -homo, Earray(k));
                }

                for (int i = 0; i < ncorrorbs; i++) {
                    if (gf_param == 2) {
                        my_gf_proj_dia(i, k - (mpi_rank * my_points)) = tmp_proj_GF(i, i);
                    } else if (gf_param == 1) {
                        my_gf_proj_dia(i, k - (mpi_rank * my_points)) = tmp_direct_GF(xbf1 - 1 + i, xbf1 - 1 + i);
                    }

                }

                if (gf_param == 2) {
                    temp_direct_Lehmann = calc_delta(tmp_proj_GF, ncorrorbs, xbf1, \
                        Over, -homo, Earray(k), 1);
                } else if (gf_param == 1) {
                    temp_direct_Lehmann = calc_delta(tmp_direct_GF, ncorrorbs, xbf1, \
                        Over_ortho_dia, -homo, Earray(k), 0);
                }
                my_D_lehmann.block(0, k - (mpi_rank * my_points), ncorrorbs, 1) = temp_direct_Lehmann;

            }
        }
        //gather data from all mpi processes
        try {
            MPI_Gather(my_D_lehmann.data(), ncorrorbs * my_points, MPI::DOUBLE_COMPLEX, D_lehmann.data(), \
                    ncorrorbs * my_points, MPI::DOUBLE_COMPLEX, 0, MPI::COMM_WORLD);
        }
        catch (MPI::Exception fail) {
            std::cerr << "MPI_Error: " << fail.Get_error_code() << " - " << fail.Get_error_string() << std::endl;
        }

        try {
            MPI_Gather(my_gf_proj_dia.data(), my_points * ncorrorbs, MPI::DOUBLE_COMPLEX, gf_proj_dia.data(), \
                    my_points * ncorrorbs, MPI::DOUBLE_COMPLEX, 0, MPI::COMM_WORLD);
        }
        catch (MPI::Exception fail) {
            std::cerr << "MPI_Error: " << fail.Get_error_code() << " - " << fail.Get_error_string() << std::endl;
        }

        //output, (only from mpi rank 0)
        if (mpi_rank == 0) {

            std::cout << "Delta computed" << std::endl;
            std::cout << "Impurity levels are: ";
            for (int i = 0; i < ncorrorbs; i++) {
                std::cout << E_imp_orbs(i) << " ";
            }
            std::cout << std::endl;

            auto t3 = std::chrono::high_resolution_clock::now();
            auto duration2 = std::chrono::duration_cast<std::chrono::seconds>(t3 - t1).count();
            std::cout << "time: " << duration2 << std::endl;

            std::ofstream outfile;
            outfile.open("delta.dat");
            for (int i = 0; i < ncorrorbs; i++) {
                for (int j = 0; j < points; j++) {
                    if (grid_param == 1) {
                        outfile << Earray(j).real() << "  " << (D_lehmann(i, j).real() - E_imp_orbs(i)) << "   " << \
                            D_lehmann(i, j).imag() << std::endl;
                    }
                    if (grid_param == 2) {
                        outfile << i << "    " << Earray(j).imag() << "  " << D_lehmann(i, j).real() - E_imp_orbs(i) << \
                                 "    " << D_lehmann(i, j).imag() << "   " << D_lehmann(i, j).real() - E_imp_orbs(i) << \
                                 D_lehmann(j, i).imag() << std::endl;
                    }
                }
                outfile << std::endl;
                outfile << std::endl;

            }
            outfile.close();

            outfile.open("gf_local.dat");
            for (int i = 0; i < ncorrorbs; i++) {
                for (int j = 0; j < points; j++) {
                    if (grid_param == 1) {
                        outfile << Earray(j).real() << "  " << gf_proj_dia(i, j).real() << "   " << \
                            gf_proj_dia(i, j).imag() << std::endl;
                    }
                    if (grid_param == 2) {
                        outfile << i << "    " << Earray(j).imag() << "  " << gf_proj_dia(j, i).real() \
                        << "    " << gf_proj_dia(j, i).imag() << std::endl;
                    }

                }
                outfile << std::endl;
                outfile << std::endl;

            }
            outfile.close();

        }
    }

    //Bath method, here hybridization for a specified subsystem can be be calculated
    if (gf_param == 0) {
        //not parallelized so only rank 0 works
        if (mpi_rank == 0) {

            //Diagonalize bath orbitals, bath = all orbitals except corr. subsystem
            Eigen::MatrixXd bath_diag;
            bath_diag = diag_bath_orbs(xbf1, xbf2, nbas, ncorrorbs, Ham_ortho, homo);
            Ham_ortho = bath_diag.transpose() * Ham_ortho_dia * bath_diag;

            //Store 1/w-e_v in Matrix G_sub, w = energy points, v = bath orbitals
            Eigen::MatrixXcd G_sub_inv = Eigen::MatrixXcd::Zero(points, (nbas - ncorrorbs));
            Eigen::MatrixXcd G_sub = Eigen::MatrixXcd::Zero(points, (nbas - ncorrorbs));

            for (int w = 0; w < points; w++) {
                for (int v = 0; v < (nbas - ncorrorbs); v++) {
                    if (v >= xbf1)
                        G_sub_inv(w, v) = Earray(w) - Ham_ortho(ncorrorbs + v, ncorrorbs + v) + homo;
                    else
                        G_sub_inv(w, v) = Earray(w) - Ham_ortho(v, v) + homo;
                    G_sub(w, v) = 1.0 / G_sub_inv(w, v);
                }
            }

            //Store V_vi^2*Sum_a(c_va^2) in V_ix_full, for decomposition truncate Sum_a by providing basis_array,
            //with orbitals to be included in part this will be stored in V_ix_part, v = bath orbitals,
            //i = impurity orbitals, a = (all/part of) bath orbitals
            Eigen::MatrixXd V_ix_part = Eigen::MatrixXd::Zero(ncorrorbs, (nbas - ncorrorbs));
            Eigen::MatrixXd V_ix_full = Eigen::MatrixXd::Zero(ncorrorbs, (nbas - ncorrorbs));
            Eigen::ArrayXd foo;
            double bar;

            for (int i = 0; i < ncorrorbs; ++i) {
                for (int k = 0; k < (nbas - ncorrorbs); k++) {

                    if (k < xbf1 - 1) {
                        if (coup_switch == 1) {
                            foo = bath_diag.col(k).array() * b_arr;
                            bar = foo.matrix().squaredNorm();
                            V_ix_part(i, k) = bar * Ham_ortho(xbf1 - 1 + i, k) * Ham_ortho(xbf1 - 1 + i, k);
                        }
                        V_ix_full(i, k) = Ham_ortho(xbf1 - 1 + i, k) * Ham_ortho(xbf1 - 1 + i, k);
                    } else {
                        if (coup_switch == 1) {
                            foo = bath_diag.col(k + ncorrorbs).array() * b_arr;
                            bar = foo.matrix().squaredNorm();
                            V_ix_part(i, k) = bar * Ham_ortho(xbf1 - 1 + i, k + ncorrorbs) \
                                                * Ham_ortho(xbf1 - 1 + i, k + ncorrorbs);
                        }
                        V_ix_full(i, k) = Ham_ortho(xbf1 - 1 + i, k + ncorrorbs) \
                                            * Ham_ortho(xbf1 - 1 + i, k + ncorrorbs);
                    }
                }

            }

            //Calculate delta
            Eigen::MatrixXcd D_partial = Eigen::MatrixXcd::Zero(points, ncorrorbs);
            Eigen::MatrixXcd D_full = Eigen::MatrixXcd::Zero(points, ncorrorbs);

            if (coup_switch == 1) {
                D_partial = G_sub * V_ix_part.transpose();
            }

            D_full = G_sub * V_ix_full.transpose();

            //Output
            std::ofstream outfile;
            if (coup_switch == 1) {
                outfile.open("delta_partial.dat");
                for (int i = 0; i < ncorrorbs; ++i) {
                    for (int k = 0; k < points; ++k) {
                        outfile << Earray(k).real() << "  " << D_partial(k, i).real() << "   " << \
                           D_partial(k, i).imag() << std::endl;
                    }
                    outfile << std::endl;
                    outfile << std::endl;
                }
                outfile.close();
            }

            outfile.open("delta_full.dat");
            for (int i = 0; i < ncorrorbs; ++i) {
                for (int k = 0; k < points; ++k) {
                    outfile << Earray(k).real() << "  " << D_full(k, i).real() << "   " << \
                           D_full(k, i).imag() << std::endl;
                }
                outfile << std::endl;
                outfile << std::endl;
            }
            outfile.close();

        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "time: " << duration << std::endl;


    MPI::Finalize();

    return (0);
}

