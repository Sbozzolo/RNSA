#!/usr/bin/env python3

# Compute Takami Rezzolla's (2011) eigenfrequencies

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 1.0
# First Stable: 05/31/17
# Last Edit: 05/31/17

import numpy as np

# Constants
a = np.zeros(6)
b = np.zeros(6)
a[0] =  2.110e-7
a[1] =  1.172e-1
a[2] = -9.599e1
a[3] =  3.621e4
a[4] = -7.757e6
a[5] =  6.978e8
b[0] =  3.357e-4
b[1] = -1.896
b[2] =  2.545e3
b[3] = -1.612e6
b[4] =  4.862e8
b[5] = -5.599e10

rhoc = float(input('Central density? '))
beta = float(input('T/W ratio? '))

terms = [a[i]*np.power(rhoc, i) + beta*b[i]*np.power(rhoc,i) for i in range(6)]
terms = [x/((4.89e-3)**2) for x in terms]
print(terms)
print("Frequency squared: {:.3E}".format(np.sum(terms)))
