#!/usr/bin/env python3

# Utility to convert polytropic units in cgs

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 1.0
# First Stable: 05/05/17
# Last Edit: 05/05/17

import numpy as np

# Constants in CGS
c     = 2.998e+10
G     = 6.674e-8
M_sun = 1.988e+33

print("DON'T USE ME! I AM BUGGY!!!")

print ("Convert polytropic units to CGS or CGS to polytropic units?")
print ("(1): Polytropic to CGS")
print ("(2): CGS to polytropic")

while True:
    user_input = input('(1) or (2)? ')

    if user_input in ['1', '2']:
        break
    else:
        print('That is not a valid option!')

if user_input == '1':
    poly_index = float(input('Polytropic index (N)? '))
    kappa = float(input('Polytropic kappa? '))
    print("What do you want to convert?")
    print("(1): Energy density")
    print("(2): Mass")
    while True:
        user_input2 = input('(1) or (2)? ')

        if user_input2 in ['1', '2']:
            break
        else:
            print('That is not a valid option!')

    if user_input2 == '1':
        # print("With K = 1 or with K = K?")
        # print("(1): K = 1 (CST)")
        # print("(2): K = K")
        # while True:
        #     user_input3 = input('(1) or (2)? ')

        #     if user_input3 in ['1', '2']:
        #         break
        #     else:
        #         print('That is not a valid option!')

        # if user_input3 == '1':
        #     pass
        # elif user_input3 == '2':
        #     pass

        energy = float(input('Value of energy (in polytropic units)? '))
        result = np.power(c, 2*poly_index)*energy/np.power(kappa, poly_index)
        print("Energy in polytropic units: ", result)
        # print("Energy in polytropic units (K = 1): ", result)
        print("Energy in CGS units: ", result)

    elif user_input2 == '2':
        mass = float(input('Value of mass (in polytropic units)? '))
        result = np.power(c/np.sqrt(G), 3)*np.power(kappa/(c*c), poly_index/2)*mass
        print("Mass in CGS units:", result)
        print("Mass in solar masses:", result/M_sun)


elif user_input == '2':
    poly_index = float(input('Polytropic index (N)? '))
    kappa = float(input('Polytropic kappa? '))
    print("What do you want to convert?")
    print("(1): Energy density")
    print("(2): Mass")
    while True:
        user_input2 = input('(1) or (2)? ')

        if user_input2 in ['1', '2']:
            break
        else:
            print('That is not a valid option!')

    if user_input2 == '1':
        energy = float(input('Value of energy (in CGS units)? '))
        result = np.power(kappa/(c*c), poly_index)*energy
        print("Energy in polytropic units: ", result)

    elif user_input2 == '2':
        massraw = float(input('Value of mass (in CGS units or solar masses)? '))
        if (massraw < 10e+20):
            mass = massraw*M_sun
        result = mass/(np.power(c/np.sqrt(G), 3)*np.power(kappa/(c*c), poly_index/2))
        print("Mass in polytropic units:", result)
