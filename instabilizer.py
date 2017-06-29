#!/usr/bin/env python3

# Instabilizer for RNS models
# Computes models on the instability line

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 2.5
# First Stable: 18/03/17
# Last Edit: 29/06/17

import argparse
import sys
import os
import shutil
import subprocess
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

def energy_cgs_to_poly(energy, index, kappa):
    return str(float(energy + "e15")/(CU_to_dens_CGS)*np.power(kappa, index))

def mass_dimensionless_to_poly(mass, index, kappa):
    return str(float(mass)/np.power(kappa, index/2))

def jmoment_dimensionless_to_poly(jmoment, index, kappa):
    return str(float(jmoment)/np.power(kappa, index))


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type = str, required = True,
                    help = "basefiles, eg. 2017_3_13_21_39")
group0 = parser.add_mutually_exclusive_group()
group0.add_argument("-p", "--poly", action="store_true",
                    help = "Use polytropic EOS")
parser.add_argument("-q", "--index", action="store", type = float,
                    help = "Polytropic index")
parser.add_argument("-k", "--kappa", action="store", type = float,
                    help = "Polytropic kappa")
group0.add_argument("-e", "--eos", type = str,
                    help = "set eos, eg. eosC")
parser.add_argument("-n", "--name", help = "name of the turning file (if not energyvsgmass)",
                    type = str)
parser.add_argument('--version', action='version', version='%(prog)s 2.5')

if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Parse arguments
eos = args.eos
poly  = args.poly
kappa = args.kappa
index = args.index

exec_path = "/home/sbozzolo/master_thesis/rns_poly/rns"


if (not eos == None):
    # Use tabulated eos
    tab_eos = True
    # Should be set!
    eos_path  = "/home/sbozzolo/master_thesis/rns4.0/EOS/" + eos
    try:
        my_file = open(eos_path)
    except IOError:
        print("EOS NOT FOUND!")
        sys.exit(2)
else:
    tab_eos = False # Polytropic
    eos = "p" + str(index)
    if (poly and (kappa == None)):
        print("You should provide the value of K!")
        sys.exit(3)

if (args.name == None):
    name = "energyvsgmass"
else:
    name = args.name

# Working dir
homedir = os.getcwd()

# Now basedir is the folder with the recap in the name
newbasedirs  = [f for f in os.listdir(os.path.join(homedir, args.file))
                if os.path.isdir(os.path.join(homedir, args.file, f))]

# Monitor the progress counting how many directories are left
j = 0

for d in newbasedirs:

    # Monitor the progress
    j +=1
    # i counts the sequences
    i = 1

    basefile =  os.path.join(homedir, args.file, d, "{}.turning".format(name))
    # Number of models
    tot = sum(1 for line in open(basefile))
    # The file is formatted in this way
    # energy gmass J/M value_of_J_or_M "A1" value_of_A_1 "A2" value_of_A_2 "B" value_of_B
    # J/M depending if the turning points is from a J or M fixed sequence
    # Stripping this data is possible to compute the model with RNS

    with open(basefile) as f:
        # Clean up previous results, which begin with the letter A
        files = [f for f in os.listdir(os.path.join(homedir, args.file)) if f.startswith("A")
                 and not os.path.isdir(os.path.join(homedir, args.file, f))]
        for ff in files: os.remove(os.path.join(homedir, args.file, ff))

        # Read file line by line
        for line in f:
            run = line.strip().split(' ')
            # RNS -s 1 means constant M, -s 2 means constant J
            # Prepare command
            if (run[2] == "STATIC"):
                if (not poly): cmd = exec_path + " -c 0.5 -d 0 -f " + eos_path + " -n 1 -r 1 -e " + run[0] + "e15"
                else:
                    cmd = exec_path + " -c 0.5 -d 0 -q poly -N " + str(index) + " -n 1 -r 1 -e " + energy_cgs_to_poly(run[0], index, kappa)
            else:
                if (run[2] == "M"):
                    task = "1"
                    if (poly): pvalue = mass_dimensionless_to_poly(run[3], index, kappa)
                if (run[2] == "J"):
                    task = "2"
                    if (poly): pvalue = jmoment_dimensionless_to_poly(run[3], index, kappa)

                if (not poly):
                    cmd = exec_path + " -c 0.5 -d 0 -f " + eos_path + " -n 1 -s " + task + " -e " + run[0] + "e15 -p " + run[3]
                else:
                    cmd = exec_path + " -c 0.5 -d 0 -q poly -N " + str(index) + " -n 1 -s " + task + " -e " +  energy_cgs_to_poly(run[0], index, kappa) + "e15 -p " + pvalue
                # Add differential rotation one parameter in case
                if (not run[5] == "0.0"): cmd += " -R diff -A {}".format(run[5])
                # Add differential rotation three parameters in case
                if (not run[7] == "0.0"): cmd += " -D {} -b {}".format(run[7], run[9])
            print (cmd)

            # Save results the try statement is required in case RNS's output value is not 0
            try:
                result = subprocess.check_output(cmd, shell = True)
                result = result.decode("utf-8")
            except Exception:
                pass

            # Progess
            print("Models: [{}/{}], Sequences [{}/{}]".format(i, tot, j, len(newbasedirs)))
            print(result)
            i += 1
            splitted = result.split()

            # If there is no result skip writing
            if(len(splitted) > 3):

                if (run[2] == "STATIC"):
                    A1value = 0.0
                    A2value = 0.0
                    Bvalue  = 0.0
                else:
                    A1value = run[5]
                    A2value = run[7]
                    Bvalue  = run[9]

                prefix_file = "inst_{}_A1_{}_A2_{}_B_{}".format(run[2], A1value, A2value, Bvalue)

                # Write result on a generic file A_value
                with open(os.path.join(homedir, args.file, d, prefix_file), 'a') as o:
                    print(result, file = o)
