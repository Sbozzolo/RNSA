#!/usr/bin/env python3

# Instabilizer for RNS models

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 1.0
# First Stable: 18/03/17
# Last Edit: 27/03/17

import argparse
import sys
import os
import shutil
import subprocess
from subprocess import Popen, PIPE

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type = str,
                    help = "basefiles, eg. 2017_3_13_21_39")
parser.add_argument("-e", "--eos", help = "set eos, eg. eosC", type = str)
parser.add_argument("-R", "--diff", help = "enable differential rotation", action = "store_true")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Parse arguments
eos = args.eos
eos_path  = "/home/sbozzolo/master_thesis/rns4.0/EOS/" + eos

try:
    my_file = open(eos_path)
except IOError:
    print("EOS NOT FOUND!")
    sys.exit(2)


homedir = "/home/sbozzolo/master_thesis/rns4.0"
basefile =  os.path.join(homedir, args.file, "energyvsgmass.turning")
# The file is formatted in this way
# energy gmass value_of_A J/M
# J/M depending if the turning points is from a J or M fixed sequence
# Stripping this data is possible to compute the model with RNS

with open(basefile) as f:
    # Clean up previous results, which begin with the letter A
    files = [f for f in os.listdir(os.path.join(homedir, args.file)) if f.startswith("A")
             and not os.path.isdir(os.path.join(homedir, args.file, f))]
    for ff in files: os.remove(os.path.join(homedir, args.file, ff))
    # Read file line by line
    for line in f:
        run = line.split(' ')
        # RNS -s 1 means constant M, -s 2 means constant J
        if (run[2] == "M"): task = "1"
        if (run[2] == "J"): task = "2"
        # Prepare command
        cmd = "./rns -d 0 -c 0.5  -f " + eos_path + " -n 1 -s " + task + " -e " + run[0] + "e15 -p " + run[3]
        # Add differential rotation in case
        if (args.diff): cmd += " -R diff -A " + run[5][0:-1]
        print (cmd)

        # Save results
        result = subprocess.check_output(cmd, shell = True)
        result = result.decode("utf-8")

        print(result)
        splitted = result.split()

        # If there is no result skip writing
        if(len(splitted) > 3):

            # Write result on a generic file A_value
            with open(os.path.join(homedir, args.file, "A" + run[5][0:-1] + run[2]), 'a') as o:
                print(result, file = o)

            # Write rmass on the same file with name appended _rmass
            with open(os.path.join(homedir, args.file, "A" + run[5][0:-1] + run[2] + "_rmass"), 'a') as o:
                print(splitted[2], file = o)

            # Write jmoment on the same file with name appended _jmoment
            with open(os.path.join(homedir, args.file, "A" + run[5][0:-1] + run[2] + "_jmoment"), 'a') as o:
                print(splitted[4], file = o)
