#!/usr/bin/env python3

# Plotter for RNS models

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 2.0
# First Stable: 13/03/17
# Last Edit: 14/04/17

import argparse
import sys
import os
import re
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import warnings
warnings.simplefilter('ignore', np.RankWarning)

cho =  ["energy", "maxenergy", "gmass", "rmass", "jmoment", "twratio",
                 "comega", "maxomega", "eomega", "radius", "rratio"]
# This is the output format
#  e_c e_max Mass Mass_0 J T/W Omega_c Omega_max Omega_e Omega_K R_e r_e grv2 grv3 r_ratio_40pcrho_c r_ratio
# dict is a dictionary to translate labels to the corresponding column number
dict = {
    "energy"      : 0,
    "maxenergy"   : 1,
    "gmass"       : 2,
    "rmass"       : 3,
    "jmoment"     : 4,
    "twratio"     : 5,
    "comega"      : 6,
    "maxomega"    : 7,
    "eomega"      : 8,
    "radius"      : 10,
    "rratio"      : 15
}
# The actual labels for the plots
labels = {
    'energy'    : "Central Energy ${\epsilon_c}\slash{c^2}~[\SI{e15}{\g\per\cm\cubed}]$",
    'maxenergy' : "Max Energy ${\epsilon_{max}}\slash{c^2}~[\SI{e15}{\g\per\cm\cubed]$",
    'gmass'     : "Gravitational Mass $M~[M_\odot]$",
    'rmass'     : "Rest Mass $M_0~[M_\odot]$",
    'jmoment'   : "Angular Momentum $\\frac{cJ}{G{M_\odot}^2}$",
    'twratio'   : "Ratio $\\frac{T}{W}$ ",
    'comega'    : "Central Angular Velocity $\Omega_c~[\si{\s^{-1}}]$",
    'maxomega'  : "Max Angular Velocity $\Omega_{max}~[\si{\s^{-1}}]$",
    'eomega'    : "Equatorial Angular Velocity $\Omega_e~[\si{\s^{-1}}]$",
    'radius'    : "Equatorial Radius $R_e~[\si{\km}]$",
    'rratio'    : "Polar Equatorial Ratio",
}

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type = str, required = True,
                    help = "basefiles, eg. \"2017_3_13_21_39\"",
                    nargs = '+')
parser.add_argument("-s", "--save", help = "save fig", action = "store_true")
parser.add_argument("-l", "--latex", help = "export to tikz", action = "store_true")
parser.add_argument("-S", "--show", help = "show fig", action = "store_true")
parser.add_argument("-r", "--rescale", help = "rescale ", action = "store_true")
parser.add_argument("-o", "--output", help = "set output folder", type = str)
parser.add_argument("-n", "--name", help = "save name", type = str)
parser.add_argument("-t", "--title", help = "set title", type = str)
parser.add_argument("-x", "--xaxis", help = "x axis, eg. energies",
                    required = True, choices = cho, type = str)
parser.add_argument("-y", "--yaxis", help = "y axis, eg. gmass",
                    required = True, choices = cho, type = str)
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()


# Set up name
if (args.name == None):
    name = '{}vs{}'.format(args.xaxis, args.yaxis)
else:
    name = args.name

# Matplotlib options, taken from
# https://github.com/MaxNoe/python-plotting/blob/master/source/siunitx_ticks.py

if plt.rcParams["text.usetex"] is False:
    plt.rcParams["text.usetex"] = True

if plt.rcParams["text.latex.unicode"] is False:
    plt.rcParams["text.latex.unicode"] = True

if "siunitx" not in plt.rcParams["text.latex.preamble"]:
    plt.rcParams["text.latex.preamble"].append(r"\usepackage{siunitx}")


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
        xlabels = [r"$\num[locale={}]{{{:.2f}}}$".format(locale, tick) for tick in xticks]
        ax.set_xticklabels(xlabels)

    if yaxis is True:
        yticks = ax.get_yticks()
        ylabels = [r"$\num[locale={}]{{{:.2f}}}$".format(locale, tick) for tick in yticks]
        ax.set_yticklabels(ylabels)

# Useful routine

# Function to find which is the best degree to fit data
def best_poly_fit(exp_x, exp_y, max_deg = 12):
    """
    Find the best degree to fit data
    Return best_degree, error and R2adj
    """
    chisq_dof = [999999]
    for deg in range(1, max_deg + 1):
        po, res, _, _, _ = np.polyfit(exp_x, exp_y, deg, full = True)
        p = np.poly1d(po)
        sum_squared = 0
        for i in range(len(exp_x)):
            sum_squared += (exp_y[i] - p(exp_x[i]))**2
        if (len(res) == 0): break;
        chisq_dof.append(sum_squared / (len(exp_x) - 2))

    best_deg = chisq_dof.index(min(chisq_dof))

    # Goodness of the fit
    RSS = sum_squared
    TSS = sum((exp_y - np.mean(exp_y))**2)
    R2 = 1 - RSS/TSS
    R2adj = R2 - (1 - R2)*(best_deg + 1)/(len(exp_x) - best_deg - 2)

    return best_deg, np.sqrt(min(chisq_dof)), R2adj


# The program actually starts here

# Working directory
homedir = os.getcwd()

# Prepare the plot
fig, ax = plt.subplots(1,1)

# Prepare colors
# The number 100 should be changed depening on how many plots are there
color = iter(cm.rainbow(np.linspace(0, 1, 100)))

# plp is an array used for the legend
plp = []

# Axis
xax = args.xaxis
yax = args.yaxis

# X and Y are arrays that contain every single point
# Used to do a global fit
X = []
Y = []

for d in args.files:
    # Base directory with current d
    basedir =  os.path.join(homedir, d)

    # Now basedir is the folder with the recap in the name
    newbasedirs  = [f for f in os.listdir(basedir)
                    if os.path.isdir(os.path.join(basedir,f))]
    oldbasedir = basedir

    for dirs in newbasedirs:

        run = d.split('_')
        eos = run[0][3:]

        files = [f for f in os.listdir(os.path.join(basedir, dirs))
                 if f.startswith("inst")]
        files.sort()

        # STATIC is the last one (if they are sorted) usually
        for f in reversed(files):
            print(d, dirs, f)
            if (f.startswith("inst_STATIC")):
                fullpath = os.path.join(basedir, dirs, f)
                static = np.genfromtxt(fullpath, max_rows = 1)
                tovx, tovy = static[dict[xax]], static[dict[yax]]
                break

        for f in files:
            fullpath = os.path.join(basedir, dirs, f)

            if (os.stat(fullpath).st_size != 0):

                if (f.startswith("inst_STATIC")):
                    pass
                else:
                    # os.system("cat {} | grep -v T= | grep -v W= >> /tmp/tmp".format(fullpath))
                    # os.system("mv /tmp/tmp {}".format(fullpath))
                    x = np.loadtxt(fullpath, usecols = dict[xax])
                    y = np.loadtxt(fullpath, usecols = dict[yax])
                    x = np.insert(x, 0, tovx)
                    y = np.insert(y, 0, tovy)
                    print(x,y)

                    if (args.rescale):
                        # EDIT HERE ACCORDING TO HOW YOU WANT TO RESCALE
                        y = y/tovy
                        x = x/tovx

                    X = np.append(X,x)
                    Y = np.append(Y,y)

                    c = next(color)
                    ax.plot(x, y, 'o', color = c)
                    run = f.split('_')
                    tmp, = plt.plot(x, y, color = c, label = "A1 = {}, A = {}, B = {}, {} const".format(run[3], run[5], run[7], run[1]))
                    plp.append(tmp)

    # Many possible fits

    # def f3(x, a, b):
    #     return 1 + a * (x**2) + b * (x**4)

    # p5 = curve_fit(f3, X, Y)
    # # print(p5)
    # xp = np.linspace(np.amin(X), np.amax(X), 100)
    # yp = [f3(x, p5[0][0], p5[0][1]) for x in xp]
    # tmp, = plt.plot(xp, yp, color = 'black', linewidth = 3.,
    #                 label = "$ {M_0}(J)= 1 + 0.51 J^2 - 0.29 J^4$")
    # plp.append(tmp)


    # poly  = np.polyfit(X, Y, 1)
    # p = np.poly1d(poly)
    # # print(poly)
    # xp = np.linspace(np.amin(X), np.amax(X), 100)
    # tmp, = plt.plot(xp, p(xp), color = 'black', linewidth = 3.,
    #                 # label = "A = {}".format(f[1:4]))
    #                             label = "$ M = 0.94 M_0 + 0.06$")
    # plp.append(tmp)

    # mmean = np.mean(m)
    # qmean = np.mean(q)
    # print("m = {} +- {}".format(mmean, np.std(m)))
    # print("q = {} +- {}".format(qmean, np.std(q)))
    # a1mean = np.mean(a1)
    # a2mean = np.mean(a2)
    # a3mean = np.mean(a3)

    # p2 = np.poly1d([a1mean,a2mean,a3mean])
    # ax.plot(xp, p2(xp), ':')
    # tmp, = plt.plot(xp, p2(xp), color = 'blue', linewidth = 3.,
    #                 # label = "A = {}".format(f[1:4]))
    #                             label = "Average quadratic")
    # plp.append(tmp)

    # def f(x, a, b):
    #     return b * np.arccosh(a * x)

    # p3 = curve_fit(f, X, Y)
    # yp = [f(x, p3[0][0], p3[0][1]) for x in xp]
    # tmp, = plt.plot(xp, yp, color = 'green', linewidth = 3.,
    #                 label = "y = b * arccosh (a x) ")
    # plp.append(tmp)

    # def f2(x, a):
    #     return a*(x-1)**(0.5)

    # p4 = curve_fit(f2, X, Y)
    # yp = [f2(x, p4[0][0]) for x in xp]
    # tmp, = plt.plot(xp, yp, color = 'purple', linewidth = 3.,
    #                 label = "y = a sqrt(x-1) ")
    # plp.append(tmp)

    # def f3(x, a, b):
    #     return 1 + a * (x**2) + b * (x**4)

    # p5 = curve_fit(f3, X, Y)
    # print(p5)
    # xp = np.linspace(np.amin(X), np.amax(X), 100)
    # yp = [f3(x, p5[0][0], p5[0][1]) for x in xp]
    # tmp, = plt.plot(xp, yp, color = 'black', linewidth = 3.,
    #                 label = "$ {M_0}(J)= 1 + 0.51 J^2 - 0.29 J^4$")
    # plp.append(tmp)

    # def f4(x, a):
    #     return  np.cosh(a * x)

    # p6 = curve_fit(f4, Y, X)
    # yp = np.linspace(np.amin(Y), np.amax(Y), 100)
    # xp = [f4(y, p6[0][0]) for y in yp]
    # tmp, = plt.plot(xp, yp, color = 'red', linewidth = 3.,
    #                 label = "x = cosh(a y) ")
    # plp.append(tmp)


    # xp = np.linspace(1,1.2, 100)
    # yp = [mmean*x + qmean for x in xp]
    # yp = [a1mean*x*x + a2mean*x + a3mean for x in xp]
    # yp = p(xp)
    # sigma = [y*0.01 for y in yp]
    # ypp = [x + y for x, y in zip(yp, sigma)]
    # ypm = [x - y for x, y in zip(yp, sigma)]
    # ypp2 = [x + 2*y for x, y in zip(yp, sigma)]
    # ypm2 = [x - 2*y for x, y in zip(yp, sigma)]
    # ypp3 = [x + 3*y for x, y in zip(yp, sigma)]
    # ypm3 = [x - 3*y for x, y in zip(yp, sigma)]

    # # tmp, = plt.plot(xp, yp, color = 'black', linewidth = 4.,
    # #                             # label = "A = {}".format(f[1:4]))
    # #                                         label = "Mean")
    # ax.fill_between(xp, ypp, ypm, alpha=.25, facecolor = 'black')
    # ax.fill_between(xp, ypp2, ypm2, alpha=.15, facecolor = 'black')
    # ax.fill_between(xp, ypp3, ypm3, alpha=.05, facecolor = 'black')
    # # plp.append(tmp)

    # tmp, = plt.plot([1, 1], [0, 0], color = 'black', linewidth = 5.,
    #                 alpha = .25, label = "\SI{1}{\percent} error")
    # plp.append(tmp)
    # tmp, = plt.plot([1, 1], [0, 0], color = 'black', linewidth = 5.,
    #                 alpha = .15, label = "\SI{2}{\percent} error")
    # plp.append(tmp)
    # tmp, = plt.plot([1, 1], [0, 0], color = 'black', linewidth = 5.,
    #                 alpha = .05, label = "\SI{3}{\percent} error")
    # plp.append(tmp)

# What to insert in legend
siunitx_ticklabels(ax)
ax.legend(handles=(plp), loc = 4)
if (not args.title == None): ax.set_title(args.title)
ax.set_ylabel(labels[args.yaxis])
ax.set_xlabel(labels[args.xaxis])

# Save to pdf
if (args.save):
    fig.tight_layout()
    if (args.output == None):
        if (len(args.files) > 1):
            fig.savefig(os.path.join("{}.pdf".format(name)), dpi = 1200)
        else:
            fig.savefig(os.path.join(basedir, "{}.pdf".format(name)), dpi = 1200)
    else:
        fig.savefig(os.path.join(args.output, "{}.pdf".format(name)), dpi = 1200)

# Save to .tikz
if (args.latex):
    from matplotlib2tikz import save as tikz_save
    if (args.output == None):
        if (len(args.files) > 1):
            tikz_save(os.path.join("{}.tikz".format(name)),
                      figureheight = '\\figureheight',
                      figurewidth = '\\figurewidth')
        else:
            tikz_save(os.path.join(basedir, "{}.tikz".format(name)),
                      figureheight = '\\figureheight',
                      figurewidth = '\\figurewidth')
    else:
        tikz_save(os.path.join(args.output, "{}.tikz".format(name)),
                  figureheight = '\\figureheight',
                  figurewidth = '\\figurewidth')

# Show via GUI
if (args.show): plt.show()
