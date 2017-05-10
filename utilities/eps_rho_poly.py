#!/usr/bin/env python3

# Utility to convert polytropic units in cgs

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 1.0
# First Stable: 05/05/17
# Last Edit: 05/05/17

import numpy as np
from scipy.optimize import newton

# Constants in CGS
c     = 2.998e+10
G     = 6.674e-8
M_sun = 1.988e+30

print ("Convert epsilon to rho or rho to epsilon?")
print ("(1): Epsilon to Rho")
print ("(2): Rho to Epslion")

while True:
    user_input = input('(1) or (2)? ')

    if user_input in ['1', '2']:
        break
    else:
        print('That is not a valid option!')

if user_input == '1':
    print("NOT SAFE!!!")
    poly_index = float(input('Polytropic index (N)? '))
    kappa = float(input('Polytropic kappa? '))
    gamma = 1 + 1/poly_index
    eps = float(input('Epsilon? '))
    def F_par(N, K, epsi):
        return lambda x : x + N*K*np.power(x, 1 + 1/N) - epsi

    F = F_par(poly_index, kappa, eps)

    print(newton(F, eps*1.01))

elif user_input == '2':
    poly_index = float(input('Polytropic index (N)? '))
    kappa = float(input('Polytropic kappa? '))
    gamma = 1 + 1/poly_index
    rho = float(input('Rho? '))
    eps = rho + poly_index*kappa*np.power(rho, gamma)
    print("Epsilon: ", eps)
