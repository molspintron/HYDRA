# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:15:48 2021

@author: Marc Philipp Bahlke (University of Hamburg)
"""

import numpy as np
import cclib as cc


'''
This is a tool which uses cclib to extract HYDRA required quantities.

The Fock matrix is evaluated from Roothaan-Hall equations, by first calculating
the overlap matrix from the following relation:
    
    C^dag*SC = 1

The Fock matrix can be calculated by using the relation:
    
    FC=SCE
    
C and E are obtained by cclib.


So far, this tool has been tested for Gaussian (closed-shell) and Turbomole 
(open and closed shell). However, it shoulkd give a general idea of how
to use other qc codes for HYDRA using cclib.

Input
------

    blanc space separated list of files for cclib
    
    For Turbomole (closed shell), e.g, <dscf.out mos coef>

    For Turbomole (open shell), e.g, <dscf.out alpha beta>

    For Gaussian, e.g,  <gaussian.log>


This tools prints
-----------------
    overlap
    hamiltonian
    eigenvalues
    coef    

    

ToDo
-----
Test open-shell case with Gaussian.
Test other QC codes than Turbomole and Gaussian.
'''


# Input
outfiles = (input('output files (blanc space separated)? '))

# Read QC output
if len(outfiles.split()) == 0:
    print('No Input provided')
elif len(outfiles.split()) == 1:
    calculation = cc.io.ccread(outfiles)
elif len(outfiles.split()) > 1:    
    outfiles = outfiles.split()
    calculation = cc.io.ccread(outfiles)

# Some empty lists
al_moenergies = []
al_mocoeffs   = []
be_moenergies = []
be_mocoeffs   = []

# Begin with reading properties
try:
    moenergies = np.asarray(calculation.moenergies)
    
    # care for number os spin components adn bring mo energies in matrix form
    # 
    if len(moenergies) == 2:
        al_moenergies = moenergies[0]
        be_moenergies = moenergies[1]
        
        siz = len(al_moenergies)
        al_zeros = np.zeros((siz,siz))
        be_zeros = np.zeros((siz,siz))
        
        for i, ene in enumerate(al_moenergies):
            al_zeros[i,i] = ene/27.2114    # !!eV to hartree!!
        
        al_moenergies = al_zeros

        for i, ene in enumerate(be_moenergies):
            be_zeros[i,i] = ene/27.2114    # !!eV to hartree!!
        
        be_moenergies = be_zeros


    else:
        al_moenergies = moenergies[0]

        siz = len(al_moenergies)
        al_zeros = np.zeros((siz,siz))

        for i, ene in enumerate(al_moenergies):
            al_zeros[i,i] = ene/27.2114    # !!eV to hartree!!
        
        al_moenergies = al_zeros
       
except:
    print('No MO energy data found')

# Get MO coefs and calculate overlap matrx according to:
# S^TCS = 1 -> C = 1/(CC^T)
# Afterwards calculate Fock matrix
try:
    mocoeffs = np.asarray(calculation.mocoeffs)
    
    if len(mocoeffs) > 1:
        
        al_mocoeffs = mocoeffs[0]
        be_mocoeffs = mocoeffs[1]
        
        inv_al_mocoeffs = np.linalg.inv(al_mocoeffs)
        inv_be_mocoeffs = np.linalg.inv(be_mocoeffs)
    elif len(mocoeffs) == 1:

        al_mocoeffs = mocoeffs[0]
        inv_al_mocoeffs = np.linalg.inv(al_mocoeffs)

    # alpha and beta overlap should be the same, so we just evaluate it once        
    overlap = np.matmul(inv_al_mocoeffs,np.matrix.transpose(inv_al_mocoeffs))
    np.savetxt('overlap', overlap, fmt='%.9e')
     
except:
    print('No MO coef data found')

# If everything is present, calculate Fock matrix F for alpha an beta
# F = SCEC^1-
if type(al_moenergies) == np.ndarray and type(al_mocoeffs) == np.ndarray:
       
    tmp1 = np.matmul(overlap, np.matrix.transpose(al_mocoeffs))# SC
    tmp2 = np.matmul(al_moenergies, np.matrix.transpose(inv_al_mocoeffs))# EC^-1

    # get Fock matrix
    F_al = np.matmul(tmp1,tmp2)
    np.savetxt('hamiltonian.1', F_al, fmt='%.9e')
    np.savetxt('eigenvalues.1', al_moenergies, fmt='%.9e')
    np.savetxt('coef.1', al_mocoeffs, fmt='%.9e')    
    
if type(be_moenergies) == np.ndarray and type(be_mocoeffs) == np.ndarray:
       
    tmp1 = np.matmul(overlap, np.matrix.transpose(be_mocoeffs))# SC
    tmp2 = np.matmul(be_moenergies, np.matrix.transpose(inv_be_mocoeffs))# EC^-1

    # get Fock matrix
    F_be = np.matmul(tmp1,tmp2)    
    np.savetxt('hamiltonian.2', F_be, fmt='%.9e')
    np.savetxt('eigenvalues.2', be_moenergies, fmt='%.9e')
    np.savetxt('coef.2', be_mocoeffs, fmt='%.9e')



    
    
