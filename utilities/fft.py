#!/usr/bin/env python3

# coding: utf-8

import numpy as np
from numpy import loadtxt as loadtxt
import matplotlib.pyplot as plt
import scipy.fftpack
import sys

try:
    data = loadtxt(sys.argv[1])
except:
    print("File not found!")
    sys.exit(1)

# Number of samplepoints
N = len(data[:,0])
x = data[:,0]
y = data[:,1]
# sample spacing
T = x[1] - x[0]

yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

plt.plot(xf, 2.0/N * np.abs(yf[0:int(N/2)]))
plt.show()
