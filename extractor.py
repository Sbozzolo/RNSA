#!/usr/bin/env python3

# Extract, to extract data from models of RNS v4

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 2.5
# First Stable: 13/03/17
# Last Edit: 15/06/17

import argparse
import sys
import os
import shutil
import re
import numpy as np

#constants, in SI
G     = 6.673e-11   # m^3/(kg s^2)
c     = 299792458   # m/s
M_sun = 1.98892e30  # kg
# Conversion factors
CU_to_km   = M_sun*G/(1000*c*c)                  # km
CU_to_s    = (M_sun*G/(c*c*c))                   # s
CU_to_dens = c*c*c*c*c*c / (G*G*G * M_sun*M_sun) # kg/m^3
CU_to_dens_CGS = CU_to_dens *1000/100**3         # g/cm^3

def energy_poly_to_cgs(energy, index, kappa, poly):
    if (not poly):
        return float(energy)*1e-15
    else:
        return float(energy)/np.power(kappa, index)*CU_to_dens_CGS*1e-15

def mass_poly_to_dimensionless(mass, index, kappa, poly):
    if (not poly):
        return float(mass)
    else:
        return float(mass)*np.power(kappa, index/2)

def jmoment_poly_to_dimensionless(jmoment, index, kappa, poly):
    if (not poly):
        return float(jmoment)
    else:
        return float(jmoment)*np.power(kappa, index)

def omega_poly_to_cgs(omega, index, kappa, poly):
    if (not poly):
        return float(omega)
    else:
        return float(omega)/np.power(kappa, index/2)*CU_to_s

def radius_poly_to_si(radius, index, kappa, poly):
    if (not poly):
        return float(radius)
    else:
        return float(radius)*np.power(kappa, index/2)*CU_to_km

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type = str, nargs = '+',
                    required = True, help = "basefiles, eg. 2017_03_13_21_39")
parser.add_argument("-a", "--all", help = "every possible data",
                    action = "store_true")
parser.add_argument("-c", "--clean", help = "clean all data files",
                    action = "store_true")
parser.add_argument("-e", "--energy", help = "central energy",
                    action = "store_true")
parser.add_argument("-em", "--maxenergy", help = "max energy",
                    action = "store_true")
parser.add_argument("-m", "--gmass", help = "gravitational mass",
                    action = "store_true")
parser.add_argument("-m0", "--rmass", help = "rest mass",
                    action = "store_true")
parser.add_argument("-j", "--jmoment", help = "angular momentum",
                    action = "store_true")
parser.add_argument("-tw", "--twratio", help = "T/W ratio",
                    action = "store_true")
parser.add_argument("-co", "--comega", help = "central angular velocity",
                    action = "store_true")
parser.add_argument("-mo", "--maxomega", help = "max angular velocity",
                    action = "store_true")
parser.add_argument("-eo", "--eomega", help = "equatorial angular velocity",
                    action = "store_true")
parser.add_argument("-R", "--radius", help = "equatorial radius",
                    action = "store_true")
parser.add_argument("-r", "--rratio", help = "polar equatorial ratio",
                    action = "store_true")
parser.add_argument("-k", "--kappa", help = "polytropic kappa",
                    action = "store", type = float)

parser.add_argument('--version', action='version', version='%(prog)s 2.0')

if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

homedir = os.getcwd()

# It can happen that there are spurious results
def is_not_ok(file):
    """Check if file contains the string 'Maximum'
    or 'sqrt"""
    with open(file, 'r') as f:
        line = f.readline()
        if line.find("Maximum") == 0:
            return True
        if line.find("sqrt") == 0:
            return True
    return False

# file counter (for stats)
i = 0

# Set alias for extract all
if (args.all):
    args.energy    = True
    args.maxenergy = True
    args.gmass     = True
    args.rmass     = True
    args.jmoment   = True
    args.twratio   = True
    args.comega    = True
    args.eomega    = True
    args.radius    = True
    args.rratio    = True

# Possible datafiles
datafiles = ["energy", "maxenergy", "gmass", "rmass", "jmoment", "twratio",
                 "comega", "maxomega", "eomega", "radius", "rratio"]

for fff in args.files:

    # Folder with only one folder inside, the one with the recap,
    # like 2017_4_10_10_14
    basedir = os.path.join(homedir, fff)

    # Now basedir is the folder with the recap in the name
    newbasedir  = [f for f in os.listdir(basedir) if os.path.isdir(os.path.join(basedir,f))][0]
    basedir = os.path.join(basedir, newbasedir)

    if ((re.match('p[0-9]\.[0-9]', newbasedir)) is not None):
        poly = True
        index = float(newbasedir[1:4])
        if (not args.clean):
            print("Found Polytrope with index", index)
            if (args.kappa == None):
                kappa = float(input('Value of K? '))
            else:
                kappa = args.kappa
    else:
        poly = False

    if ("static" in basedir):
        static = True
    else:
        static = False

    if (args.clean):
        # If it is static the structure is simple
        if (static):
            # Iterate over possible data files, defined in datafiles array
            for df in datafiles:
                if (os.path.exists(os.path.join(basedir, df + ".dat"))):
                    os.remove(os.path.join(basedir, df + ".dat"))
                    print("Removed " + os.path.join(basedir, df + ".dat"))
                    i += 1
        # If it is not static the structure is complex
        else:
            # folders with A1 fixed
            foldersA1  = [f for f in os.listdir(basedir) if os.path.isdir(os.path.join(basedir, f))]
            foldersA1.sort()
            for fA1 in foldersA1:
                # folders with A2 fixed
                foldersA2  = [f for f in os.listdir(os.path.join(basedir, fA1)) if os.path.isdir(os.path.join(basedir, fA1, f))]
                foldersA2.sort()
                for fA2 in foldersA2:
                    # folders with beta fixed
                    foldersB  = [f for f in os.listdir(os.path.join(basedir, fA1, fA2)) if os.path.isdir(os.path.join(basedir, fA1, fA2, f))]
                    foldersB.sort()
                    for fB in foldersB:
                        # folders with different values of jmoment or rmass
                        folders  = [f for f in os.listdir(os.path.join(basedir, fA1, fA2, fB))  if os.path.isdir(os.path.join(basedir, fA1, fA2, fB, f))]
                        folders.sort()
                        # inside these folders there are the actual data
                        # Iterate over folders
                        for d in folders:
                            # Iterate over data files
                            for df in datafiles:
                                if (os.path.exists(os.path.join(basedir, fA1, fA2, fB, d, df + ".dat"))):
                                    os.remove(os.path.join(basedir, fA1, fA2, fB, d, df + ".dat"))
                                    print("Removed " + os.path.join(basedir, fA1, fA2, fB, d, df + ".dat"))
                                    i += 1
    # If it is not clean
    else:
        if (static):
            files  = [f for f in os.listdir(basedir) if f.endswith(".out") and
                      os.stat(os.path.join(basedir, f)).st_size != 0 and
                      sum(1 for line in open(os.path.join(basedir, f))) < 4 and
                      not is_not_ok(os.path.join(basedir, f))]
            files.sort()
            if (args.energy):
                energypath = os.path.join(basedir, "energy.dat")
                energyfile = open(energypath, "w")
            if (args.maxenergy):
                maxenergypath = os.path.join(basedir, "maxenergy.dat")
                maxenergyfile = open(maxenergypath, "w")
            if (args.gmass):
                gmasspath = os.path.join(basedir, "gmass.dat")
                gmassfile = open(gmasspath, "w")
            if (args.rmass):
                rmasspath = os.path.join(basedir, "rmass.dat")
                rmassfile = open(rmasspath, "w")
            if (args.jmoment):
                jmomentpath = os.path.join(basedir, "jmoment.dat")
                jmomentfile = open(jmomentpath, "w")
            if (args.twratio):
                twratiopath = os.path.join(basedir, "twratio.dat")
                twratiofile = open(twratiopath, "w")
            if (args.comega):
                comegapath = os.path.join(basedir, "comega.dat")
                comegafile = open(comegapath, "w")
            if (args.maxomega):
                maxomegapath = os.path.join(basedir, "maxomega.dat")
                maxomegafile = open(maxomegapath, "w")
            if (args.eomega):
                eomegapath = os.path.join(basedir, "eomega.dat")
                eomegafile = open(eomegapath, "w")
            if (args.radius):
                radiuspath = os.path.join(basedir, "radius.dat")
                radiusfile = open(radiuspath, "w")
            if (args.rratio):
                rratiopath = os.path.join(basedir, "rratio.dat")
                rratiofile = open(rratiopath, "w")
            for f in files:
                    with open(os.path.join(basedir, f), 'r') as ff:
                        data = ff.readline().strip().split()
                        if (args.energy):    print(energy_poly_to_cgs(data[0], index, kappa, poly),
                                                   file = energyfile)
                        if (args.maxenergy): print(energy_poly_to_cgs(data[1], index, kappa, poly),
                                                   file = maxenergyfile)
                        if (args.gmass):     print(mass_poly_to_dimensionless(data[2], index, kappa, poly),
                                                   file = gmassfile)
                        if (args.rmass):     print(mass_poly_to_dimensionless(data[3], index, kappa, poly),
                                                   file = rmassfile)
                        if (args.jmoment):   print(jmoment_poly_to_dimensionless(data[4], index, kappa, poly),
                                                   file = jmomentfile)
                        if (args.twratio):   print(data[5],  file = twratiofile)
                        if (args.comega):    print(omega_poly_to_cgs(data[6], index, kappa, poly),
                                                   file = comegafile)
                        if (args.maxomega):  print(omega_poly_to_cgs(data[7], index, kappa, poly),
                                                   file = maxomegafile)
                        if (args.eomega):    print(omega_poly_to_cgs(data[8], index, kappa, poly),
                                                   file = eomegafile)
                        if (args.radius):    print(radius_poly_to_si(data[10], index, kappa, poly),
                                                   file = radiusfile)
                        if (args.rratio):    print(data[15], file = rratiofile)
            if (args.energy):
                energyfile.close()
                print("Written ", energypath)
                i += 1
            if (args.maxenergy):
                maxenergyfile.close()
                print("Written ", maxenergypath)
                i += 1
            if (args.gmass):
                gmassfile.close()
                print("Written ", gmasspath)
                i += 1
            if (args.rmass):
                rmassfile.close()
                print("Written ", rmasspath)
                i += 1
            if (args.jmoment):
                jmomentfile.close()
                print("Written ", jmomentpath)
                i += 1
            if (args.twratio):
                twratiofile.close()
                print("Written ", twratiopath)
                i += 1
            if (args.comega):
                comegafile.close()
                print("Written ", comegapath)
                i += 1
            if (args.maxomega):
                maxomegafile.close()
                print("Written ", maxomegapath)
                i += 1
            if (args.eomega):
                eomegafile.close()
                print("Written ", eomegapath)
                i += 1
            if (args.radius):
                radiusfile.close()
                print("Written ", radiuspath)
                i += 1
            if (args.rratio):
                rratiofile.close()
                print("Written ", rratiopath)
                i += 1
        # If it is not static the structure is complex
        else:
            # folders with A1 fixed
            foldersA1  = [f for f in os.listdir(basedir) if os.path.isdir(os.path.join(basedir, f))]
            foldersA1.sort()
            for fA1 in foldersA1:
                # folders with A2 fixed
                foldersA2  = [f for f in os.listdir(os.path.join(basedir, fA1))  if os.path.isdir(os.path.join(basedir, fA1, f))]
                foldersA2.sort()
                for fA2 in foldersA2:
                    # folders with beta fixed
                    foldersB  = [f for f in os.listdir(os.path.join(basedir, fA1, fA2)) if os.path.isdir(os.path.join(basedir, fA1, fA2, f))]
                    foldersB.sort()
                    for fB in foldersB:
                        # folders with different values of jmoment or rmass
                        path_B = os.path.join(basedir, fA1, fA2, fB)
                        folders  = [f for f in os.listdir(path_B)  if os.path.isdir(os.path.join(basedir, fA1, fA2, fB, f))]
                        folders.sort()
                        # inside these folders there are the actual data
                        # Iterate over folders
                        for d in folders:
                            # files is the array which contains the datafiles
                            path_folder = os.path.join(path_B, d)
                            files  = [f for f in os.listdir(path_folder) if f.endswith(".out") and
                                      os.stat(os.path.join(path_folder, f)).st_size != 0 and
                                      sum(1 for line in open(os.path.join(path_folder, f))) < 4 and
                                      not is_not_ok(os.path.join(path_folder, f))]
                            files.sort()
                            if (args.energy):
                                energypath = os.path.join(path_folder, "energy.dat")
                                energyfile = open(energypath, "w")
                            if (args.maxenergy):
                                maxenergypath = os.path.join(path_folder, "maxenergy.dat")
                                maxenergyfile = open(maxenergypath, "w")
                            if (args.gmass):
                                gmasspath = os.path.join(path_folder, "gmass.dat")
                                gmassfile = open(gmasspath, "w")
                            if (args.rmass):
                                rmasspath = os.path.join(path_folder, "rmass.dat")
                                rmassfile = open(rmasspath, "w")
                            if (args.jmoment):
                                jmomentpath = os.path.join(path_folder, "jmoment.dat")
                                jmomentfile = open(jmomentpath, "w")
                            if (args.twratio):
                                twratiopath = os.path.join(path_folder, "twratio.dat")
                                twratiofile = open(twratiopath, "w")
                            if (args.comega):
                                comegapath = os.path.join(path_folder, "comega.dat")
                                comegafile = open(comegapath, "w")
                            if (args.maxomega):
                                maxomegapath = os.path.join(path_folder, "maxomega.dat")
                                maxomegafile = open(maxomegapath, "w")
                            if (args.eomega):
                                eomegapath = os.path.join(path_folder, "eomega.dat")
                                eomegafile = open(eomegapath, "w")
                            if (args.radius):
                                radiuspath = os.path.join(path_folder, "radius.dat")
                                radiusfile = open(radiuspath, "w")
                            if (args.rratio):
                                rratiopath = os.path.join(path_folder, "rratio.dat")
                                rratiofile = open(rratiopath, "w")
                            for f in files:
                                with open(os.path.join(path_folder, f), 'r') as ff:
                                    data = ff.readline().strip().split()
                                    if (args.energy):    print(energy_poly_to_cgs(data[0], index, kappa, poly),
                                                               file = energyfile)
                                    if (args.maxenergy): print(energy_poly_to_cgs(data[1], index, kappa, poly),
                                                               file = maxenergyfile)
                                    if (args.gmass):     print(mass_poly_to_dimensionless(data[2], index, kappa, poly),
                                                               file = gmassfile)
                                    if (args.rmass):     print(mass_poly_to_dimensionless(data[3], index, kappa, poly),
                                                               file = rmassfile)
                                    if (args.jmoment):   print(jmoment_poly_to_dimensionless(data[4], index, kappa, poly),
                                                               file = jmomentfile)
                                    if (args.twratio):   print(data[5],  file = twratiofile)
                                    if (args.comega):    print(omega_poly_to_cgs(data[6], index, kappa, poly),
                                                               file = comegafile)
                                    if (args.maxomega):  print(omega_poly_to_cgs(data[7], index, kappa, poly),
                                                               file = maxomegafile)
                                    if (args.eomega):    print(omega_poly_to_cgs(data[8], index, kappa, poly),
                                                               file = eomegafile)
                                    if (args.radius):    print(radius_poly_to_si(data[10], index, kappa, poly),
                                                               file = radiusfile)
                                    if (args.rratio):    print(data[15], file = rratiofile)
                            if (args.energy):
                                energyfile.close()
                                print("Written ", energypath)
                                i += 1
                            if (args.maxenergy):
                                maxenergyfile.close()
                                print("Written ", maxenergypath)
                                i += 1
                            if (args.gmass):
                                gmassfile.close()
                                print("Written ", gmasspath)
                                i += 1
                            if (args.rmass):
                                rmassfile.close()
                                print("Written ", rmasspath)
                                i += 1
                            if (args.jmoment):
                                jmomentfile.close()
                                print("Written ", jmomentpath)
                                i += 1
                            if (args.twratio):
                                twratiofile.close()
                                print("Written ", twratiopath)
                                i += 1
                            if (args.comega):
                                comegafile.close()
                                print("Written ", comegapath)
                                i += 1
                            if (args.maxomega):
                                maxomegafile.close()
                                print("Written ", maxomegapath)
                                i += 1
                            if (args.eomega):
                                eomegafile.close()
                                print("Written ", eomegapath)
                                i += 1
                            if (args.radius):
                                radiusfile.close()
                                print("Written ", radiuspath)
                                i += 1
                            if (args.rratio):
                                rratiofile.close()
                                print("Written ", rratiopath)
                                i += 1
if (args.clean):
    print("Done! Removed {} files".format(i))
else:
    print("Done! Written {} files".format(i))


#  e_c e_max Mass Mass_0 J T/W Omega_c Omega_max Omega_e Omega_K R_e r_e grv2 grv3 r_ratio_40pcrho_c r_ratio
