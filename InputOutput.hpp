#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

Eigen::MatrixXd read_denmat(std::fstream& in, int n, int rank, std::string name);

std::vector<float> read_input(std::fstream& in, int rank);

Eigen::ArrayXd read_intvector(std::fstream& in, int n, int rank);

void print_output(int nbas, int xbf1, int xbf2, int grid_param, int beta, int points, float step, float smear, double homo, int unit, int diag_param, int gf_param, int coup_switch, int my_points, int mpi_size, int omp_size);

#endif
