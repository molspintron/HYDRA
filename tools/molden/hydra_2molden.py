# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:15:48 2021

Written by Dr. Marc Philipp Bahlke
University of Hamburg

This is a small script which reads the header of a 
molden file and uses HYDRA output files to create
impurity and bath molden files.

Requires:
--------

molden.input     (molden file of a converged QC calculation)
molden.bath      (HYDRA output)
molden.impurity  (HYDRA output)

Output:
------
molden_bath.input
molden_impurity.input

"""


# First read molden.input, molden.bath and molden.impurity

molden_header =  []
switch = True

try:
    with open('molden.input') as file:
        for line in file:
            if switch == True:
                molden_header.append(line.split())
            if '[MO]' in line:
                switch = False
except OSError:
    print('Can not find molden.input. Please use tm2molden.')
    
bath = []

try:
    with open('molden.bath') as file:
        for line in file:
            bath.append(line.split())

except OSError:
    print('Can not find molden.bath')
    
    
impurity = []

try:
    with open('molden.impurity') as file:
        for line in file:
            impurity.append(line.split())

except OSError:
    print('Can not find molden.impurity')

    
# Now print molden readable files

outdat = 'molden_bath.input'
    
with open(outdat, 'w') as file:
    for line in molden_header:
        if any(x == '[GTO]' for x in line):
            file.write('[5D]\n') # required for real valued spherical harmonics           
        file.write(''.join( x+' ' for x in line )+'\n')
    for line in bath:
        file.write(''.join( x+' ' for x in line )+'\n')

outdat = 'molden_impurity.input'

with open(outdat, 'w') as file:
    for line in molden_header:
        if any(x == '[GTO]' for x in line):
            file.write('[5D]\n') # required for real valued spherical harmonics
        file.write(''.join( x+' ' for x in line )+'\n')
    for line in impurity:
        file.write(''.join( x+' ' for x in line )+'\n')




