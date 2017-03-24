#!/usr/bin/env python3

# Plotter for RNS models

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 1.5
# First Stable: 13/03/17
# Last Edit: 23/03/17

import argparse
import sys
import os
import re
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib2tikz import save as tikz_save
import warnings
warnings.simplefilter('ignore', np.RankWarning)

# from framework import siunitx_ticklabels
from framework import best_poly_fit

# Parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type = str, required = True,
                    help = "basefiles, eg. 2017_3_13_21_39")
parser.add_argument("-x", "--xaxis", help = "x axis, eg. energies",
                    required = True, type = str)
parser.add_argument("-y", "--yaxis", help = "y axis, eg. gmass",
                    required = True, type = str)
parser.add_argument("--static", help = "set static folder",
                    type = str)
parser.add_argument("-s", "--save", help = "save fig", action = "store_true")
parser.add_argument("-S", "--show", help = "show fig", action = "store_true")
parser.add_argument("--nokepler", help = "disable kepler solution", action = "store_true")
group1 = parser.add_mutually_exclusive_group()
group1.add_argument("--noturning", help = "disable turning points", action = "store_true")
group1.add_argument("-o", "--output", help = "produce a numerical output", action = "store_true")

parser.add_argument("-l", "--latex", help = "export to PGF", action = "store_true")
parser.add_argument("-n", "--name", help = "save name", type = str)
parser.add_argument("-t", "--title", help = "set title", type = str)

parser.add_argument('--version', action='version', version='%(prog)s 1.5')

if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Parse arguments

fig, ax = plt.subplots(1,1)

# basedir is where data are

homedir = "/home/sbozzolo/master_thesis/rns4.0"
basedir =  os.path.join(homedir, args.files)

# Set a name if is not provided

if (args.name == None):
    name = '{}vs{}'.format(args.xaxis, args.yaxis)
else:
    name = args.name

# Don't comput turning points
notur = args.noturning
if (notur == None): notur = False
# Numerical output
out = args.output
if (out):
    outpath = os.path.join(basedir, "{}.output".format(name))
    outfile = open(outpath, "w")
    tourpath = os.path.join(basedir, "{}.turning".format(name))
    tourfile = open(tourpath, "w")

# Very very very not sleek!

# If static is provided do not search for static in the folder ma in args.static
if (not args.static == None):
    folderstaticUP  = [f for f in os.listdir(args.static) if re.search(r'static',f)]
    folderstatic  = [f for f in os.listdir(os.path.join(homedir, args.static, folderstaticUP[0])) if re.search(r'static',f)]
    fullpathstaticx = os.path.join(homedir, args.static, folderstaticUP[0], folderstatic[0], args.xaxis + ".dat")
    fullpathstaticy = os.path.join(homedir, args.static, folderstaticUP[0], folderstatic[0], args.yaxis + ".dat")
else:
    folderstaticUP  = [f for f in os.listdir(basedir) if re.search(r'static',f)]
    folderstatic  = [f for f in os.listdir(os.path.join(basedir, folderstaticUP[0])) if re.search(r'static',f)]
    fullpathstaticx = os.path.join(basedir, folderstaticUP[0], folderstatic[0], args.xaxis + ".dat")
    fullpathstaticy = os.path.join(basedir, folderstaticUP[0], folderstatic[0], args.yaxis + ".dat")

xstatic = np.loadtxt(fullpathstaticx)
ystatic = np.loadtxt(fullpathstaticy)

# Critical point for the static sequence
if (not notur):
    # Polynomial fit
    polystatic = np.polyfit(xstatic, ystatic, best_poly_fit(xstatic, ystatic)[0])
    ps = np.poly1d(polystatic)
    # Find maxima
    crit = ps.deriv().r
    r_crit = crit[crit.imag == 0].real
    test = ps.deriv(2)(r_crit)
    xs_crit = r_crit[test < 0]
    ys_crit = ps(xs_crit)
    if (len(xs_crit) == 0):
        xs_crit = [0]
        ys_crit = [0]

if (not args.nokepler):
    folderkeplerUP  = [f for f in os.listdir(basedir) if re.search(r'kepler',f)]
    folderkepler  = [f for f in os.listdir(os.path.join(basedir, folderkeplerUP[0])) if re.search(r'kepler',f)]
    fullpathkeplerx = os.path.join(basedir, folderkeplerUP[0], folderkepler[0], args.xaxis + ".dat")
    fullpathkeplery = os.path.join(basedir, folderkeplerUP[0], folderkepler[0], args.yaxis + ".dat")
    xkepler = np.loadtxt(fullpathkeplerx)
    ykepler = np.loadtxt(fullpathkeplery)

# folderkepler  = [f for f in os.listdir(homedir) if f.startswith(args.files) and re.search(r'kepler',f)]
# fullpathkeplerx = os.path.join(homedir, folderkepler[0], args.xaxis + ".dat")
# fullpathkeplery = os.path.join(homedir, folderkepler[0], args.yaxis + ".dat")
# xkepler = np.loadtxt(fullpathkeplerx)
# ykepler = np.loadtxt(fullpathkeplery)

# Plot static sequence and kepler one
static, = plt.plot(xstatic, ystatic, "-", color = "black", label = 'Nonrotating sequence')
if (not args.nokepler): kepler, = plt.plot(xkepler, ykepler, "--",  color = "red", label = 'Mass shedding limit')

# Dealing with data

pltcounter = 0

# foldersA are every folder in basedir but static and kepler, with various values of A
foldersA  = [f for f in os.listdir(basedir) if re.search(r'^((?!static).)*$',f)
             and re.search(r'^((?!kepler).)*$',f) and os.path.isdir(os.path.join(basedir,f))]
# Better sort, you never know
foldersA.sort()

# Every folder has it's own critical point, prepare the array
critplt = []

if plt.rcParams["text.usetex"] is False:
    plt.rcParams["text.usetex"] = True
    # print("\nWARNING: text.usetex is now set to True\n")

if plt.rcParams["text.latex.unicode"] is False:
    plt.rcParams["text.latex.unicode"] = True
    # print("\nWARNING: text.latex.unicode is now set to True\n")

if "siunitx" not in plt.rcParams["text.latex.preamble"]:
    plt.rcParams["text.latex.preamble"].append(r"\usepackage{siunitx}")

if "microtype" not in plt.rcParams["text.latex.preamble"]:
    plt.rcParams["text.latex.preamble"].append(r"\usepackage[protrusion=true,factor=900]{microtype}")


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

for dA in foldersA:
    # One color for every curve
    color = iter(cm.rainbow(np.linspace(0, 1, len(foldersA))))

    # folders  = [f for f in os.listdir(os.path.join(basedir, dA))]
    # folders are the folders with fixed A (not static, kepler and ign) ign means ignore
    folders  = [f for f in os.listdir(os.path.join(basedir, dA)) if re.search(r'^((?!static).)*$',f)
                and re.search(r'^((?!kepler).)*$',f) and re.search(r'^((?!ign).)*$',f)]
    folders.sort()

    # From the name it is possible to extract useful info
    run = dA.split('_')

    # Add here for other tasks
    if ("rmass" in dA):   task = "rmass"
    if ("jmoment" in dA): task = "jmoment"

    if (not notur):
        critx = np.array(xs_crit)
        crity = np.array(ys_crit)

    if(out and not notur):
        print("Turning point for the static sequence")
        print("x = {}, y = {}".format(xs_crit[0], ys_crit[0]))
        print("Turning point for the static sequence", file = outfile)
        print("x = {}, y = {}".format(xs_crit[0], ys_crit[0]), file = outfile)
        # print("{} {}".format(xs_crit, ys_crit), file = tourfile)

    for d in folders:
        # d are the folders which contains actual data
        fullpathx = os.path.join(basedir, dA, d, args.xaxis + ".dat")
        fullpathy = os.path.join(basedir, dA, d, args.yaxis + ".dat")
        # Ignore if a file is empty
        if (os.stat(fullpathx).st_size != 0 and os.stat(fullpathy).st_size != 0):
            x = np.loadtxt(fullpathx)
            y = np.loadtxt(fullpathy)
            ax.plot(x, y, ',')

            # Extract parameters
            run2 = d.split('_')

            if(out and not notur):
                print("Turning point for the sequence A = {}, {} = {}".format(run[3], run[8][0], run2[7][1:5]))
                print("Turning points for the sequence A = {}, {} = {}".format(run[3], run[8][0], run2[7][1:5]), file = outfile)

            # Perform interpolation
            if (not notur):
                # Poly fit
                poly = np.polyfit(x, y, best_poly_fit(x,y)[0])
                p = np.poly1d(poly)
                crit = p.deriv().r
                r_crit = crit[crit.imag == 0].real
                test = p.deriv(2)(r_crit)
                if (task == "rmass"):
                    x_crit = r_crit[test > 0]
                if (task == "jmoment"):
                    x_crit = r_crit[test < 0]
                y_crit = p(x_crit)

                xp = np.linspace(np.amin(np.amin(x)), np.amax(x), 100)
                xx = [x for x in xp if p(x) >= ps(x)]
                # ax.plot(xx, p(xx), ':')

                # Number of possible turning points
                zeros = 0

                for k in range(len(x_crit)):
                    # Critical points should be discarded if they are out of boundaries, if they are too close to
                    # a boundary (which means that maybe they are spourious) and if they are below the static line
                    if (ps(x_crit[k]) < p(x_crit[k]) and x_crit[k] > np.amin(x) and x_crit[k] < np.amax(x)
                        and x_crit[k] < np.amax(x) - 0.03*np.amax(x) and x_crit[k] > np.amin(x) + 0.03*np.amin(x)):
                        critx = np.append(critx, x_crit[k])
                        crity = np.append(crity, y_crit[k])
                        zeros += 1
                        if (out):
                            print("x = {}, y = {}".format(x_crit[k], y_crit[k]))
                            print("x = {}, y = {}".format(x_crit[k], y_crit[k]), file = outfile)
                            print("{} {} {} {} {} {}".format(x_crit[k], y_crit[k],
                                                             run[8][0], run2[7][1:5], "A", run[3]), file =
                                  tourfile)

                if (zeros > 1): print("WARNING: MULTIPLE TURNING POINTS")

    if (out and notur):
        print("")
        print("", file = outfile)

    # print(critx, crity)

    # Fit critical points
    if (not notur):
        # Critical fit
        crit_deg, crit_chi, R2adj = best_poly_fit(critx, crity, 4)
        polycrit = np.polyfit(critx, crity, crit_deg) #, cov = True)
        if (out):
            print("A = {}, {} constant".format(run[3], run[8][0]), file = outfile)
            print("Fitted with a Polynom of degree {}".format(crit_deg), file = outfile)
            print("Residual standard error of the fit {:.2E}".format(crit_chi), file = outfile)
            print("R2 adjusted {:.4f}".format(R2adj), file = outfile)
            print("Coefficients:")
            # for k in range(len(polycrit)):
                # print("Coefficient of x^{}: {} +- {}".format(len(polycrit)- 1 - k,
                                                             # polycrit[k], np.sqrt(cov[k][k])), file = outfile)
            print("", file = outfile)
            print("Fitted with a Polynom of degree {}".format(crit_deg))
            print("Residual standard error of the fit {:.2e}".format(crit_chi))
            print("R2 adjusted {:.4f}".format(R2adj))
            print("Coefficients:")
            # for k in range(len(polycrit)):
                # print("Coefficient of x^{}: {} +- {}".format(len(polycrit)- 1 - k, polycrit[k], np.sqrt(cov[k][k])))
            print("")

        # Prepare a plot
        pc = np.poly1d(polycrit)
        xp = np.linspace(np.amin(critx), np.amax(critx), 100)

        # c is the color of the polt
        c = next(color)
        tmp, = plt.plot(xp, pc(xp), '-', color = c,
                                        label = "A = {}, {} constant".format(run[3], run[8][0]))
        ax.plot(critx, crity, 'o', color = c)

        critplt.append(tmp)

    pltcounter += 1

# Plotting the plot, these are the labels
labels = {
    'energy'    : "Central Energy ${\epsilon_c}\slash{c^2}~[\SI{e15}{\g\cm^{-3}}]$",
    'maxenergy' : "Max Energy ${\epsilon_{max}}\slash{c^2}~[\SI{e15}{\g\cm^{-3}}]$",
    'gmass'     : "Gravitational Mass $M~[M_\odot]$",
    'rmass'     : "Rest Mass $M_0~[M_\odot]$",
    'jmoment'   : "Angular Momentum $\frac{cJ}{G{M_\odot}^2}$",
    'twratio'   : "Ratio T W ",
    'comega'    : "Central Angular Velocity $\Omega_c~[\si{\s^{-1}}]$",
    'maxomega'  : "Max Angular Velocity $\Omega_{max}~[\si{\s^{-1}}]$",
    'eomega'    : "Equatorial Angular Velocity $\Omega_e~[\si{\s^{-1}}]$",
    'radius'    : "Equatorial Radius $R_e~[\si{\km}]$",
    'rratio'    : "Polar Equatorial Ratio",
}

# siunitx_ticklabels(ax)
ax.set_xlabel(labels[args.xaxis])
ax.set_ylabel(labels[args.yaxis])
ax.set_title(args.title)

# What to insert in legend
plops = [static]
if (not notur): plops = plops + critplt
if (not args.nokepler): plops = plops + [kepler]

ax.legend(handles=(plops) , loc = 4)

if (args.save):
    fig.tight_layout()
    fig.savefig(os.path.join(basedir, name + ".pdf"), format = 'pdf', dpi = 1200)

if (args.latex):
    tikz_save(os.path.join(basedir, name + ".tex"), dpi = 600)

if (out):
    print("Output produced in {}".format(outpath))
    if (not notur): print("Turning points in {}".format(tourpath))
    outfile.close()
    if (not notur): tourfile.close()

if (args.show): plt.show()
