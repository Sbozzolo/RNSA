#!/usr/bin/env python3

# Plotter for RNS models

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 2.1
# First Stable: 13/03/17
# Last Edit: 16/05/17

import argparse
import sys
import os
import re
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import warnings
# Ignore polyfit warning
warnings.simplefilter('ignore', np.RankWarning)

# Parse arguments

cho =  ["energy", "maxenergy", "gmass", "rmass", "jmoment", "twratio",
                 "comega", "maxomega", "eomega", "radius", "rratio"]

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type = str, help = "basefiles, eg. 2017_3_13_21_39")
parser.add_argument("-x", "--xaxis", help = "x axis, eg. energies",
                    required = True, choices = cho, type = str)
parser.add_argument("-y", "--yaxis", help = "y axis, eg. gmass",
                    required = True, choices = cho, type = str)
parser.add_argument("--static", help = "set static folder",
                    type = str, required = True)
parser.add_argument("-s", "--save", help = "save fig", action = "store_true")
parser.add_argument("-S", "--show", help = "show fig", action = "store_true")
group1 = parser.add_mutually_exclusive_group()
group1.add_argument("--noturning", help = "disable turning points", action = "store_true")
group1.add_argument("-o", "--output", help = "produce a numerical output", action = "store_true")
parser.add_argument("-l", "--latex", help = "export to PGF", action = "store_true")
parser.add_argument("-n", "--name", help = "save name", type = str)
parser.add_argument("-t", "--title", help = "set title", type = str)

parser.add_argument('--version', action='version', version='%(prog)s 2.0')

if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()


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
        xlabels = [r"$\num[locale={}]{{{:.2E}}}$".format(locale, tick) for tick in xticks]
        ax.set_xticklabels(xlabels)

    if yaxis is True:
        yticks = ax.get_yticks()
        ylabels = [r"$\num[locale={}]{{{:.2E}}}$".format(locale, tick) for tick in yticks]
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

homedir = os.getcwd()

fig, ax = plt.subplots(1,1)

# Set a name if is not provided

if (args.name == None):
    name = '{}vs{}'.format(args.xaxis, args.yaxis)
else:
    name = args.name

# Don't comput turning points (can be useful when plotting not energy-gmass but other
# quantities)

notur = args.noturning
if (notur == None): notur = False

# Numerical output
out = args.output

# Plot static solution
folderstatic  = [f for f in os.listdir(args.static) if re.search(r'static',f)]
fullpathstaticx = os.path.join(homedir, args.static, folderstatic[0], args.xaxis + ".dat")
fullpathstaticy = os.path.join(homedir, args.static, folderstatic[0], args.yaxis + ".dat")

xstatic = np.loadtxt(fullpathstaticx)
ystatic = np.loadtxt(fullpathstaticy)

# Critical point for the static sequence
if (not notur):
    # Polynomial fit
    polystatic = np.polyfit(xstatic, ystatic, best_poly_fit(xstatic, ystatic, 10)[0])
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

# Plot static sequence
static, = plt.plot(xstatic, ystatic, "-", color = "black", label = 'Nonrotating sequence')

if (not args.files == None):

    # basedir is where data are

    basedir = os.path.join(homedir, args.files)

    # Dealing with data

    # Every folder has it's own critical point, prepare the array
    critplt = []

    # Now basedir is the folder with the recap in the name
    newbasedirs  = [f for f in os.listdir(basedir) if os.path.isdir(os.path.join(basedir,f))]
    oldbasedir = basedir

    # Total number of expected sequences
    ntot = 0

    for dirr in newbasedirs:
        ntot += int((dirr.strip().split('_'))[-1])

    # One color for every curve
    color = iter(cm.rainbow(np.linspace(0, 1, ntot)))

    for dirs in newbasedirs:

        basedir = os.path.join(oldbasedir, dirs)

        if (out):
            # Outfile contains the location of the turning points in easily readable
            # format (is basically the STDOUT)
            outpath = os.path.join(basedir, "{}.output".format(name))
            outfile = open(outpath, "w")
            # Tourfile contains a table with detailed information about the turning points
            # To use with instabilizer
            tourpath = os.path.join(basedir, "{}.turning".format(name))
            tourfile = open(tourpath, "w")

        # Add here for other tasks
        if ("rmass" in basedir):   task = "rmass"
        if ("jmoment" in basedir): task = "jmoment"

        # folders with A1 fixed
        foldersA1  = [f for f in os.listdir(basedir)
                      if os.path.isdir(os.path.join(basedir, f))]
        foldersA1.sort()
        for fA1 in foldersA1:
            # folders with A2 fixed
            foldersA2  = [f for f in os.listdir(os.path.join(basedir, fA1))
                          if os.path.isdir(os.path.join(basedir, f, fA1))]
            foldersA2.sort()
            for fA2 in foldersA2:
                # folders with beta fixed
                foldersB  = [f for f in os.listdir(os.path.join(basedir, fA1, fA2))
                             if os.path.isdir(os.path.join(basedir, f, fA1, fA2))]
                foldersB.sort()
                for fB in foldersB:
                    # folders with different values of jmoment or rmass
                    folders  = [f for f in os.listdir(os.path.join(basedir, fA1, fA2, fB))
                                if os.path.isdir(os.path.join(basedir, f, fA1, fA2, fB))]
                    folders.sort()
                    # inside these folders there are the actual data

                    if (not notur):
                        # Prepare the array which contains the critial points with the static one
                        critx = np.array(xs_crit)
                        crity = np.array(ys_crit)

                    if(out and not notur):
                        print("Turning point for the static sequence")
                        print("x = {}, y = {}".format(xs_crit[0], ys_crit[0]))
                        print("Turning point for the static sequence", file = outfile)
                        print("x = {}, y = {}".format(xs_crit[0], ys_crit[0]), file = outfile)
                        print("{} {} {}".format(xs_crit[0], ys_crit[0], "STATIC"), file = tourfile)


                    # Iterate over folders
                    base_path = os.path.join(basedir, fA1, fA2, fB)
                    for d in folders:
                        # d are the folders which contains actual data
                        fullpathx = os.path.join(base_path, d, args.xaxis + ".dat")
                        fullpathy = os.path.join(base_path, d, args.yaxis + ".dat")


                        # Ignore if a file is empty
                        if (os.stat(fullpathx).st_size != 0 and os.stat(fullpathy).st_size != 0
                            and sum(1 for line in open(fullpathx)) > 1 and  sum(1 for line in open(fullpathy)) > 1):
                            x = np.loadtxt(fullpathx)
                            y = np.loadtxt(fullpathy)
                            ax.plot(x, y, ',')

                            # Extract parameters
                            run = d.split('_')

                            if(out and not notur):
                                print("A1 = {}, A2 = {}, B = {}, {} = {}".format(run[2], run[4], run[6], run[11], run[12]))
                                print("A1 = {}, A2 = {}, B = {}, {} = {}".format(run[2], run[4], run[6], run[11], run[12]),
                                      file = outfile)

                            # Perform interpolation
                            if (not notur):
                                # Poly fit
                                poly = np.polyfit(x, y, best_poly_fit(x,y)[0], 15)
                                p = np.poly1d(poly)
                                crit = p.deriv().r
                                r_crit = crit[crit.imag == 0].real
                                test = p.deriv(2)(r_crit)
                                if (task == "rmass"):
                                    # Minima
                                    x_crit = r_crit[test > 0]
                                if (task == "jmoment"):
                                    # Maxima
                                    x_crit = r_crit[test < 0]
                                y_crit = p(x_crit)

                                # Interpolated points
                                xp = np.linspace(np.amin(np.amin(x)), np.amax(x), 100)
                                xx = [x for x in xp if p(x) >= ps(x)]
                                # ax.plot(xx, p(xx), ':')

                                # Number of possible turning points
                                zeros = 0

                                for k in range(len(x_crit)):
                                    # Critical points should be discarded if they are out of boundaries, if they are too close to
                                    # a boundary (which means that maybe they are spourious) and if they are below the static line
                                    if (ps(x_crit[k]) < p(x_crit[k]) and x_crit[k] > np.amin(x) and x_crit[k] < np.amax(x)
                                        and x_crit[k] < np.amax(x) - 0.03*np.amax(x) and x_crit[k] > np.amin(x) + 0.03*np.amin(x)
                                        and sum(1 for line in open(fullpathx)) > 25 and sum(1 for line in open(fullpathy)) > 25
                                        and (np.abs(ps(x_crit[k]) - p(x_crit[0])) > 0.05*ys_crit[0] or np.abs(x_crit[k] - xs_crit[0]) < 0.1*xs_crit[0])
                                    ):
                                        zeros += 1

                                if(zeros > 1):
                                    print("WARNING: MULTIPLE TURNING POINTS")

                                    # If multiple turning points are found let the user decide
                                    while True:
                                        user_input = input('{}, {} (1) or {} {} (2) or discard (3) '.format(x_crit[0],
                                                                                                            y_crit[0], x_crit[1], y_crit[1]))
                                        if user_input in ['1', '2', '3']:
                                            break
                                        else:
                                            print('That is not a valid option!')

                                    if user_input == '1':
                                        critx = np.append(critx, x_crit[0])
                                        crity = np.append(crity, y_crit[0])
                                        if (out):
                                            print("x = {}, y = {}".format(x_crit[0], y_crit[0]))
                                            print("x = {}, y = {}".format(x_crit[0], y_crit[0]), file = outfile)
                                            print("{} {} {} {} {} {} {} {} {} {}".format(x_crit[0], y_crit[0],
                                                                                      run[11], run[12], "A1",
                                                                                      run[2], "A2", run[4],
                                                                                      "B", run[6]), file = tourfile)
                                    elif user_input == '2':
                                        critx = np.append(critx, x_crit[1])
                                        crity = np.append(crity, y_crit[1])
                                        if (out):
                                            print("x = {}, y = {}".format(x_crit[1], y_crit[1]))
                                            print("x = {}, y = {}".format(x_crit[1], y_crit[1]), file = outfile)
                                            print("{} {} {} {} {} {} {} {} {} {}".format(x_crit[1], y_crit[1],
                                                                                      run[11], run[12], "A1",
                                                                                      run[2], "A2", run[4],
                                                                                      "B", run[6]), file = tourfile)
                                    # Ignore the point (safest solution)
                                    elif user_input == '3':
                                        pass
                                # Single zero
                                else:
                                    # I don't know which one it is
                                    for k in range(len(x_crit)):
                                    # Critical points should be discarded if they are out of boundaries, if they are too close to
                                    # a boundary (which means that maybe they are spourious) and if they are below the static line
                                        if (ps(x_crit[k]) < p(x_crit[k]) and x_crit[k] > np.amin(x) and x_crit[k] < np.amax(x)
                                            and x_crit[k] < np.amax(x) - 0.03*np.amax(x) and x_crit[k] > np.amin(x) + 0.03*np.amin(x)
                                            and sum(1 for line in open(fullpathx)) > 25 and sum(1 for line in open(fullpathy)) > 25
                                            and (np.abs(ps(x_crit[k]) - p(x_crit[0])) > 0.05*ys_crit[0] or np.abs(x_crit[k] - xs_crit[0]) < 0.1*xs_crit[0])
                                            ):
                                            critx = np.append(critx, x_crit[k])
                                            crity = np.append(crity, y_crit[k])
                                            if (out):
                                                print("x = {}, y = {}".format(x_crit[k], y_crit[k]))
                                                print("x = {}, y = {}".format(x_crit[k], y_crit[k]), file = outfile)
                                                print("{} {} {} {} {} {} {} {} {} {}".format(x_crit[k], y_crit[k],
                                                                                      run[11], run[12], "A1",
                                                                                      run[2], "A2", run[4],
                                                                                      "B", run[6]), file = tourfile)

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
                            print("A1 = {}, A2 = {}, B = {}, {} constant".format(run[2], run[4], run[6], run[11]), file = outfile)
                            print("Fitted with a Polynom of degree {}".format(crit_deg), file = outfile)
                            print("Residual standard error of the fit {:.2E}".format(crit_chi), file = outfile)
                            print("R2 adjusted {:.4f}".format(R2adj), file = outfile)
                            # print("Coefficients:")
                            # for k in range(len(polycrit)):
                                # print("Coefficient of x^{}: {} +- {}".format(len(polycrit)- 1 - k,
                                                                             # polycrit[k], np.sqrt(cov[k][k])), file = outfile)
                            print("", file = outfile)
                            print("Fitted with a Polynom of degree {}".format(crit_deg))
                            print("Residual standard error of the fit {:.2e}".format(crit_chi))
                            print("R2 adjusted {:.4f}".format(R2adj))
                            # print("Coefficients:")
                            # for k in range(len(polycrit)):
                                # print("Coefficient of x^{}: {} +- {}".format(len(polycrit)- 1 - k, polycrit[k], np.sqrt(cov[k][k])))
                            print("")

                        # Prepare a plot
                        pc = np.poly1d(polycrit)
                        xp = np.linspace(np.amin(critx), np.amax(critx), 100)

                        # c is the color of the polt
                        c = next(color)
                        ax.plot(xp, pc(xp), '-', color = c)

                        if (task == "jmoment"):
                            tmp, = plt.plot(critx, crity, 'o', color = c,
                                            label = "{}, A1 = {}, A2 = {}, B = {}, {} constant".format(run[0], run[2], run[4],
                                                                                                       run[6], run[11]))
                        elif (task == "rmass"):
                            tmp, = plt.plot(critx, crity, '*', color = c,
                                            label = "{}, A1 = {}, A2 = {}, B = {}, {} constant".format(run[0], run[2], run[4],
                                                                                                       run[6], run[11]))

                        critplt.append(tmp)


# Plotting the plot, these are the labels
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


G     = 6.673e-11   # m^3/(kg s^2)
c     = 299792458   # m/s
M_sun = 1.98892e30  # kg

CU_to_dens = c*c*c*c*c*c / (G*G*G * M_sun*M_sun) # kg/m^3
CU_to_dens_CGS = CU_to_dens *1000/100**3         # g/cm^3

ax2 = ax.twiny()

if (args.xaxis == "energy"):
    new_tick_locations = ax.get_xticks()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels((new_tick_locations)/CU_to_dens_CGS*1e18)
    ax2.set_xlabel("Dimensionless Central Energy $[\num{1e-3}]$")

if (args.xaxis == "maxenergy"):
    new_tick_locations = ax.get_xticks()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels((new_tick_locations)/CU_to_dens_CGS*1e18)
    ax2.set_xlabel("Dimensionless Max Energy $[\num{1e-3}]$")

if (args.yaxis == "energy"):
    new_tick_locations = ax.get_yticks()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(new_tick_locations)
    ax2.set_yticklabels((new_tick_locations)/CU_to_dens_CGS*1e18)
    ax2.set_ylabel("Dimensionless Central Energy $[\num{1e-3}]$")

if (args.xaxis == "maxenergy"):
    new_tick_locations = ax.get_yticks()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(new_tick_locations)
    ax2.set_yticklabels((new_tick_locations)/CU_to_dens_CGS*1e18)
    ax2.set_ylabel("Dimensionless Max Energy $[\num{1e-3}]$")

siunitx_ticklabels(ax)
ax.set_xlabel(labels[args.xaxis])
ax.set_ylabel(labels[args.yaxis])
if (not args.title == None):
    ax.set_title(args.title)

# What to insert in legend
plops = [static]
if (not notur): plops = plops + critplt

ax.legend(handles=(plops) , loc = 4)

# Save the pdf
if (args.save):
    fig.tight_layout()
    fig.savefig(os.path.join(basedir, name + ".pdf"), format = 'pdf', dpi = 1200)

# Save the tikz file
if (args.latex):
    from matplotlib2tikz import save as tikz_save
    tikz_save(os.path.join(basedir, name + ".tikz"),
              figureheight = '\\figureheight',
              figurewidth = '\\figurewidth')

# Save output and turning
if (out):
    print("Output produced in {}".format(outpath))
    if (not notur): print("Turning points in {}".format(tourpath))
    outfile.close()
    if (not notur): tourfile.close()

# Show with GUI
if (args.show): plt.show()
