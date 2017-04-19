#!/usr/bin/env python3

# Utility to have a TIKZ rainbow

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 1.0
# First Stable: 19/04/17
# Last Edit: 19/04/17

import argparse
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib2tikz import save as tikz_save

parser = argparse.ArgumentParser(description='Produce a TIKZ rainbow')
parser.add_argument('points', metavar='N', type=int,
                    help='how many colors do you need?')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

color = iter(cm.rainbow(np.linspace(0, 1, args.points)))

fig, ax = plt.subplots(1,1)

for i in range(args.points):
    c = next(color)
    ax.plot([0],[0], color = c)

    # Trash everything to /dev/null
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")

    try:
        tikz_save("/tmp/tmp.tex")
    finally:
        # Restore stdout
        sys.stdout.close()
        sys.stdout = old_stdout

# Read and grep definecolo
with open('/tmp/tmp.tex', 'r') as f:
    for line in f:
        line = line.rstrip()
        if "definecolor" in line:
            print (line)

# Cleanup
os.remove("/tmp/tmp.tex")
