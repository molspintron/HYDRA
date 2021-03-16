#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <sstream>
#include <complex>
#include "InputOutput.hpp"

//Reads basis_array
Eigen::ArrayXd read_intvector(std::fstream& in, int n, int rank)
{
    int value;
    std::string line;
    Eigen::ArrayXd b_arr = Eigen::ArrayXd::Zero(n);
    if(!in.is_open()){
        std::cerr << "MPI_rank " << rank << ": basis_array did not open" << std::endl;}
    while (in >> value)
    {
        b_arr(value-1) = 1.0;
    }
    return b_arr;
}


//Reads matrices
Eigen::MatrixXd read_denmat(std::fstream& in, int n, int rank, std::string name)
{
    //cout << "test1\n";
    //fstream out;
    int i, j;
    std::string buf;
    Eigen::MatrixXd mat(n,n);
    if(!in.is_open()){
        std::cerr << "MPI_rank " << rank << ": " << name << " did not open" << std::endl;}

    for(i = 0; i<n; i++){
        for(j=0; j<n; j++){
            in >> buf;
            auto fixD(buf.find_first_of("Dd"));
            if(fixD != std::string::npos){
                buf[fixD] = 'E';
            }
            mat(i,j) = std::stod(buf);
        }
    }
    in.close();
    return mat;

}

//Reads input parameters from delta.global.in
std::vector<float> read_input(std::fstream& in, int rank)
{
    std::vector<float> param;
    int i = 0;
    std::string line;
    if (!in.is_open()){
        std::cerr << "MPI_rank " << rank <<  ": delta.global.in did not open" << std::endl;}
    while (getline(in, line))
    {
        std::istringstream strm(line);
        float value;
        strm >> value;
        param.push_back(value);
        i++;
    }
    in.close();
    return param;

}

void print_output(int nbas, int xbf1, int xbf2, int grid_param, int beta, int points, float step, float smear, \
    double homo, int unit, int diag_param, int gf_param, int coup_switch, int my_points, int mpi_size, int omp_size)
{
    std::cout << "*********************************" << std::endl;
    std::cout << "*                               *" << std::endl;
    std::cout << "*               HYDRA           *" << std::endl;
    std::cout << "*          Version:2.0.0        *" << std::endl;
    std::cout << "*                               *" << std::endl;
    std::cout << "*********************************" << std::endl;

    std::cout << "Basis functions: " << nbas << std::endl;
    std::cout << "Impurity orbitals: " << xbf1 << " - " << xbf2 << std::endl;
    if(gf_param != 0) {
    std::cout << "Parallelization: " << mpi_size << " MPI_ranks with " << omp_size << " OpenMP threads each";
    std::cout << std::endl << "                 " << my_points / omp_size << "energy points per thread, " << points
              << \
              " total points" << std::endl;
    }
    std::cout << "Energy grid: " << points;
    if (grid_param == 1) {
        std::cout << "real axis, " << step << "eV steps, " << smear << "eV smear" << std::endl;
    } else if (grid_param == 2) {
        std::cout << "matsubara axis, " << beta << " beta" << std::endl;
    } else {
        std::cerr << "ERROR: Energy grid, 1: real,  2: Matsubara, input: " << \
                grid_param << std::endl;
        exit(EXIT_FAILURE);
    }

        std::cout << "HOMO: " << homo << "eV" << std::endl;
       
        std::cout << "Diagonalization: ";
        if (diag_param == 1) {
            std::cout << "on" << std::endl;
        } else if (diag_param == 0) {
            std::cout << "off" << std::endl;
        } else {
            std::cerr << "ERROR: Diagonalization, 0: off,  1: on, input: " << \
                    diag_param << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Delta unit: ";
        if (unit == 1) {
            std::cout << "hartree" << std::endl;
        } else if (unit == 2) {
            std::cout << "eV" << std::endl;
        } else {
            std::cerr << "ERROR: Delta unit, 1: hartree,  2: eV, input: " << \
                    unit << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Method for GF: ";
        if (gf_param == 0) {
            std::cout << "diagonalize bath hamiltonian (NO INVERSION, fast)" << std::endl;
        } else if (gf_param == 1) {
            std::cout << "GF = (wI-H)^(-1) (INVERSION, slow)" << std::endl;
        } else if (gf_param == 2) {
            std::cerr << "ERROR: Lehmann representation and projection to subspace does not work yet!" \
            << std::endl;
            exit(EXIT_FAILURE);
        } else {
            std::cerr << "ERROR: GF method: 0: diag. bath, 1: (wI-H)^(-1) (NOT RECOMMENDED), input" << \
                    gf_param << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Calculate hybridization of part: ";
        if (coup_switch == 1) {
            std::cout << "on" << std::endl;
        } else if (coup_switch == 0) {
            std::cout << "off" << std::endl;
        } else {
            std::cerr << "ERROR: HF^(part), 1: on,  0: off, input: " << \
                    coup_switch << std::endl;
            exit(EXIT_FAILURE);
        }

}


