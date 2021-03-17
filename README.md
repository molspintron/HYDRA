--------------------------------
HYDRA
-------------------------------

HYDRA is a tool to compute hybridization functions using atom-centered basis sets.
It can perform local decompositions of hybridization functions.

The tool was developed by Dr. Marc Philipp Bahlke (in Fortran)
and re-written in C++ by Michaela Schneeberger, for an 
improved performance.


############
Installation
############

1. Download and unzip the Eigen library (http://eigen.tuxfamily.org/)
2. Create 'build' directory (mkdir build)
3. a) For omp version:
    cp source/main_omp.cpp to ./build/main.cpp and source/Makefile_omp to ./build/Makefile
    To switch off omp, add the option -fopenmp to the compile command in the Makefile. This will result in a serial version.
3. b) Copy source/InputOutput.cpp and source/InputOutput.hpp to ./build
3. c) cd ./build
3. d) Include the path to Eigen library in the Makefile & Edit Makefile according to your systems settings
3. e) type 'make' in folder with Hydra to compile it
4. a) For mpi version:
    cd ./build
4. b) cmake -G "Unix Makefiles" ../source
4. c) make all


############
Usage
############

HYDRA requires:

-Input-file         (filename=delta.global.in)
-Overlap matrix     (filename=overlap)
-Hamiltonian matrix (filename=hamiltonian)
-Eigenvalues        (filename=eigenvalues)
-MO coefficients    (filename=coef)
-Basis set indices for local decomposition (filename=basis_array)

delta.global.in looks like the following:

------------------------------------------------------------------------------------------------------------
        3000                     ! nbas, total number of basis functions
        2948                     ! xbf1, index of 1st basis function of impurity 
        2952                     ! xbf2, index of last basis function of impurity
        1                        ! 1= real, 2= matsu; type of energy axis, matsu is Matsubara axis
        100                      ! beta, inverse temperature, only relevant for Matsubara axis
        1001                     ! points, number of data points for the hybridization function
        0.02                     ! step size
        0.1                      ! smear,  smearing or value of the imaginary off-set
        -3.959                   ! Fermi energy
        1                        ! unit (0 = hartree, 1 = eV)
        1                        ! diagonalize (true if 1), note that this is not always required, check if impurity block is already diagonal!
        0                        ! method, 0 = bath, 1 = use Lehmann, 2 = projector (not fully implemented yet)
        1                        ! 1 = calculate a decomposition (only for method 0, requiers basis_array input)
------------------------------------------------------------------------------------------------------------

For Turbomole, the required matrices can be obtained by
a tool which is provided in ./tools/turbomole. It is required
to remove the file extensions in some cases, according to the
required files mentioned above. 

E.g. eigenvalues.dat -> eigenvalues
     hamiltonian.1 -> hamiltonian
     etc.

Usage for other QC codes:
Basically, you just need the required matrices as mentioned above in the exact same format (see ./examples).



##########
Plotting environment (bath) orbitals and impurity orbitals
##########

Hydra outputs 'truncated' molden files. To make them readable by molden,
it is required to add the 'header' of a full molden file to the truncated
ones. This can be done by hand or by a python script provided in:

 ./tools/molden/hydra_2molden.py

This script needs the truncated molden files provided by hydra and in addition
a full molden file which can, e.g. for Turbomole, created with the script
tm2molden (part of Turbomole).


Notes/TODOs:

- The impurity orbital indices are not allowed to start from 1. This is, e.g., the case
  when  one of the impurity atoms is the 1st atom in the xyz/coord file and the s functions
  are included in the local impurity sub-space.


