#!/usr/bin/env python3

# Plotter for RNS models

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 2.0
# First Stable: 13/03/17
# Last Edit: 11/04/17

import argparse
import sys
import os
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import cm
import numpy as np
import warnings
# Suppress Polyfit warning
warnings.simplefilter('ignore', np.RankWarning)

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type = str, required = True,
                    help = "basefiles, eg. \"2017_3_13_21_39\"",
                    nargs = '+')
parser.add_argument("-o", "--output", help = "set output file",
                    type = str, required = True)

parser.add_argument('--version', action='version', version='%(prog)s 2.0')

if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Matplotlib options, taken from
# https://github.com/MaxNoe/python-plotting/blob/master/source/siunitx_ticks.py

if plt.rcParams["text.usetex"] is False:
    plt.rcParams["text.usetex"] = True
    # print("\nWARNING: text.usetex is now set to True\n")

if plt.rcParams["text.latex.unicode"] is False:
    plt.rcParams["text.latex.unicode"] = True
    # print("\nWARNING: text.latex.unicode is now set to True\n")

if "siunitx" not in plt.rcParams["text.latex.preamble"]:
    plt.rcParams["text.latex.preamble"].append(r"\usepackage{siunitx}")

if "sisetup" not in plt.rcParams["text.latex.preamble"]:
    plt.rcParams["text.latex.preamble"].append(r"\sisetup{inter-unit-product =  \ensuremath{{}\;{}}}")


def siunitx_ticklabels(ax=None, locale="US", xaxis=True, yaxis=True):
    """
    This function uses siunitx to create the ticklabels
    Main reason is for adjusting the decimal marker properly.
    The function takes 4 arguments:
        ax=None     the matplotlib axes to operate on
                    if set to None (Standard) this will be the current axes
        locale="DE" The locale parameter for siunitx, one of
                    "UK", "US", "DE", "FR" oder "ZA"
        xaxis=True  Boolean, if True the labels for the xaxis are set
        yaxis=True  Boolean, if True the labels for the yaxis are set
    """

    if ax is None:
        ax = plt.gca()

    if xaxis is True:
        xticks = ax.get_xticks()
        xlabels = [r"$\num[locale={}]{{{:.2E}}}$".format(locale, tick) for tick in xticks]
        ax.set_xticklabels(xlabels)

    if yaxis is True:
        yticks = ax.get_yticks()
        ylabels = [r"$\num[locale={}]{{{:.2E}}}$".format(locale, tick) for tick in yticks]
        ax.set_yticklabels(ylabels)

labels = [
    "Central Energy ${\epsilon_c}\slash{c^2}~[\SI{e15}{\g\cm^{-3}}]$",
    "Max Energy ${\epsilon_{max}}\slash{c^2}~[\SI{e15}{\g\cm^{-3}}]$",
    "Gravitational Mass $M~[M_\odot]$",
    "Rest Mass $M_0~[M_\odot]$",
    "Angular Momentum $\\frac{cJ}{G{M_\odot}^2}$",
    "Ratio \\frac{T}{W} ",
    "Central Angular Velocity $\Omega_c~[\si{\s^{-1}}]$",
    "Max Angular Velocity $\Omega_{max}~[\si{\s^{-1}}]$",
    "Equatorial Angular Velocity $\Omega_e~[\si{\s^{-1}}]$",
    "Kepler Angular Velocity $\Omega_{kep}~[\si{\s^{-1}}]$",
    "Equatorial Radius $R_e~[\si{\km}]$",
    "Polar Equatorial Ratio",
]
quantities = [
    "Central Energy",
    "Max Energy",
    "Gravitational Mass",
    "Rest Mass",
    "Angular Momentum",
    "Ratio T W ",
    "Central Angular Velocity",
    "Max Angular Velocity",
    "Equatorial Angular Velocity",
    "Kepler Angular Velocity",
    "Equatorial Radius",
    "Polar Equatorial Ratio",
]
#  e_c e_max Mass Mass_0 J T/W Omega_c Omega_max Omega_e Omega_K R_e r_e grv2 grv3 r_ratio_40pcrho_c r_ratio


# Current directory
homedir = os.getcwd()

# Every plot should stay in a different page
with PdfPages(args.output) as pdf:
    for d in args.files:
        # Now basedir is the folder with the recap in the name
        newbasedirs  = [f for f in os.listdir(os.path.join(homedir, d)) if os.path.isdir(os.path.join(homedir, d, f))]

        for dd in newbasedirs:
            # dd is the folder with the recap in the name
            basedir =  os.path.join(homedir, d, dd)

            run = dd.split('_')
            eos = run[0]

            # Loops over pairs of quantities
            for i in range(11):
                for j in range(i, 11):
                    if (not i == j):
                        fig, ax = plt.subplots(1,1)

                        # Treat differently the static file
                        staticfile = [f for f in os.listdir(basedir) if f.startswith("inst_STATIC")][0]

                        # Extract only the relevant quantities
                        # I used genfromtext because the static file right now has also debus info on
                        # the second and third line, so I use max_rows
                        tovdata =   np.genfromtxt(os.path.join(basedir, staticfile), dtype = 'float',
                                                  usecols = (0,1,2,3,4,5,6,7,8,9,10,15),
                                                  max_rows = 1)

                        files  = [f for f in os.listdir(basedir) if (f.startswith("inst_J")
                                                                     or f.startswith("inst_M"))]
                        files.sort()

                        # For the legend
                        plp = []

                        # One color for every curve
                        color = iter(cm.rainbow(np.linspace(0, 1, len(files))))
                        # For monitoring the progress
                        q = 1
                        for f in files:
                            # f is like inst_J_A1_1.0_A2_1.5_B_0.8

                            fullpathdata = os.path.join(basedir, f)
                            if (os.stat(fullpathdata).st_size != 0
                                and sum(1 for line in open(os.path.join(fullpathdata))) > 2):
                                print(f)
                                data = np.loadtxt(fullpathdata, dtype = 'float',
                                                  usecols = (0,1,2,3,4,5,6,7,8,9,10,15))

                                x = np.append(tovdata[i], data[:,i])
                                y = np.append(tovdata[j], data[:,j])
                                print("Working on {}, {}, [{}/{}]".format(quantities[i],quantities[j], q, len(files)))
                                q +=1

                                c = next(color)
                                fs = f.split('_')
                                tmp, = plt.plot(x, y, 'o',
                                                label = "EOS {}, A1 = {}, A2 = {}, B = {}, {} const".format(eos, fs[3], fs[5], fs[7], fs[1]),
                                                color = c)

                                # Fit with a poly
                                poly = np.polyfit(x, y, 4)
                                p = np.poly1d(poly)

                                xp = np.linspace(np.amin(np.amin(x)), np.amax(x), 100)
                                ax.plot(xp, p(xp), ':', color = c)

                                # tmp, = plt.plot(xp, p(xp), color = pltcolours[pltcounter], # linewidth = 3.,
                                                            # label = "EOS {}, {} constant".format(eos, f[-1]))

                                plp.append(tmp)

                        # Plot the plot
                        siunitx_ticklabels(ax)
                        ax.set_title("{} vs {}".format(quantities[j], quantities[i]))
                        ax.set_xlabel(labels[i])
                        ax.set_ylabel(labels[j])

                        # what to insert in legend
                        ax.legend(handles=(plp), loc = 4)

                        fig.tight_layout()
                        # dpi are low to speed up the process
                        pdf.savefig(fig, dpi = 150)
