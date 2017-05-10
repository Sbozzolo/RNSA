#!/usr/bin/env python3

# Utility to convert dimensionless units in cgs

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 1.0
# First Stable: 08/05/17
# Last Edit: 08/05/17

import numpy as np

# Constants in CGS
c     = 2.998e+10
G     = 6.674e-8
M_sun = 1.988e+33

print ("Convert dimensionless units to CGS or CGS to dimensionless units?")
print ("(1): Dimensionless to CGS")
print ("(2): CGS to dimensionless")

while True:
    user_input = input('(1) or (2)? ')

    if user_input in ['1', '2']:
        break
    else:
        print('That is not a valid option!')

if user_input == '1':
    print("What do you want to convert?")
    print("(1): Energy densities")
    print("(2): Lengths")
    while True:
        user_input2 = input('(1) or (2)? ')

        if user_input2 in ['1', '2']:
            break
        else:
            print('That is not a valid option!')

    if user_input2 == '1':
        energy = float(input('Value of energy (in dimensionless units)? '))
        result = np.power(c, 6)/(np.power(G, 3)*np.power(M_sun, 2))*energy
        print("Energy in dimensionless units: ", energy)
        print("Energy in g/cm^3: ", result)

    elif user_input2 == '2':
        length = float(input('Value of length (in dimensionless units)? '))
        result = G*M_sun/(c*c)*length/100000
        print("Leght in dimensionless units:", length)
        print("Length in km:", result)


elif user_input == '2':
    print("What do you want to convert?")
    print("(1): Energy densities")
    print("(2): Lengths")
    while True:
        user_input2 = input('(1) or (2)? ')

        if user_input2 in ['1', '2']:
            break
        else:
            print('That is not a valid option!')

    if user_input2 == '1':
        energy = float(input('Value of energy (in g/cm^3)? '))
        result = np.power(np.power(c, 6)/(np.power(G, 3)*np.power(M_sun, 2)), -1)*energy
        print("Energy in dimensionless units: ", result)
        print("Energy in g/cm^3: ", energy)

    elif user_input2 == '2':
        length = float(input('Value of length (in km)? '))
        result = np.power(G*M_sun/(c*c)/100000, -1)*length
        print("Leght in dimensionless units:", result)
        print("Length in km:", length)
