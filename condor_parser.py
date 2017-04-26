#!/usr/bin/env python3

# Condor Parser for RNS v4.0 and HTCondor 8.6.1
# Prepare a bunch of condor.sumbit and than submit them

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 2.0
# First Stable: 13/03/17
# Last Edit: 11/04/17

import argparse
import sys
import os
import shutil
import subprocess
from subprocess import Popen, PIPE
import time
import datetime

# Parse command line arguments

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--eos", type = str, required = True,
                    help = "set eos, eg. eosC")
parser.add_argument("-t", "--task", type = str, choices = ["static", "rmass", "jmoment"],
                    required = True, help = "define which kind of sequences have to be computed")
parser.add_argument("-o", "--out", type = str,
                    help = "set output folder")
parser.add_argument('-N', "--nsequences", action="store", type = int, help =
                    "set number of sequences")
parser.add_argument('-n', "--nmodels", action="store", type = int,
                    required = True, help = "set number of models for each sequence")
parser.add_argument('-na', "--namodels", action="store", type = int,
                    help = "set number of fixed values of A (one-parameter)")
parser.add_argument('-e1', "--energy1", action="store", type = float,
                    required = True, help = "initial central energy")
parser.add_argument('-e2', "--energy2", action="store", type = float,
                    required = True, help = "final central energy")
parser.add_argument('-A1', "--rotation1", action="store", type = float,
                    help = "set initial rotation parameter A")
parser.add_argument('-A2', "--rotation2", action="store", type = float,
                    help = "set final rotation parameter A")
group1 = parser.add_mutually_exclusive_group()
group1.add_argument("-m1", "--rmass1", action="store", type = float,
                    help = "initial rest mass (in M_sun)")
group1.add_argument("-j1", "--jmoment1", action="store", type = float,
                    help = "initial angular momentum (in GM_sun^2/c)")
group2 = parser.add_mutually_exclusive_group()
group2.add_argument("-m2", "--rmass2", action="store", type = float,
                    help = "final rest mass (in M_sun)")
group2.add_argument("-j2", "--jmoment2", action="store", type = float,
                    help = "final angular momentum (in GM_sun^2/c)")
group3 = parser.add_mutually_exclusive_group()
group3.add_argument("-R", "--diff", help = "enable differential rotation",
                    action = "store_true")
group3.add_argument("-R3", "--diff3", help = "enable differential rotation with 3 parameters",
                    action = "store_true")
parser.add_argument('-B1', "--beta1", action="store", type = float,
                    help = "set initial beta")
parser.add_argument('-B2', "--beta2", action="store", type = float,
                    help = "set final beta")
parser.add_argument('-A1i', "--A1initial", action="store", type = float,
                    help = "set initial A1")
parser.add_argument('-A1f', "--A1final", action="store", type = float,
                    help = "set final A1")
parser.add_argument('-A2i', "--A2initial", action="store", type = float,
                    help = "set initial A2")
parser.add_argument('-A2f', "--A2final", action="store", type = float,
                    help = "set final A2")
parser.add_argument('-na1', "--namodels1", action="store", type = int,
                    help = "set number of fixed values of A1")
parser.add_argument('-na2', "--namodels2", action="store", type = int,
                    help = "set number of fixed values of A2")
parser.add_argument('-nB', "--nbmodels", action="store", type = int,
                    help = "set number of fixed values of beta")
parser.add_argument("-c", "--condor", help = "sumbit jobs to Condor",
                    action = "store_true")
parser.add_argument('--version', action='version', version='%(prog)s 2.0')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

task = args.task
eos  = args.eos
n    = args.nmodels
na   = args.namodels
N    = args.nsequences
e1   = args.energy1
e2   = args.energy2
m1   = args.rmass1
m2   = args.rmass2
j1   = args.jmoment1
j2   = args.jmoment2
ai   = args.rotation1
af   = args.rotation2
b1   = args.beta1
b2   = args.beta2
nb   = args.nbmodels
na1  = args.namodels1
na2  = args.namodels2
a1initial  = args.A1initial
a2initial  = args.A2initial
a1final  = args.A1final
a2final  = args.A2final

# End parsing CLI arguments

exec_path = "/home/sbozzolo/master_thesis/rns4.0/rns"
eos_path  = "/home/sbozzolo/master_thesis/rns4.0/EOS/" + eos

try:
    my_file = open(eos_path)
except IOError:
    print("EOS NOT FOUND!")
    sys.exit(2)

# Define a new folder based on time if -o is not provided

# basefolder is where data will be stored, full path

if (args.out != None):
    basefolder = args.out
else:
    now = datetime.datetime.now()
    basefolder = "{}_{}_{}_{}_{}".format(now.year, now.month, now.day, now.hour, now.minute)

# Add absolute path
basefolder = os.path.join(os.getcwd(), basefolder)

# Remove old results, if they exist
if os.path.exists(basefolder):
    shutil.rmtree(basefolder)
os.mkdir(basefolder)

# Energy step = final energy - initial energy / number of models
step = (e2 - e1)/(n - 1)

# Treat separatly the static case
if (task == "static"):
    parentfolder = os.path.join(basefolder, "{}_Ei_{:.2e}_Ef_{:.2e}_{}/".format(eos, e1, e2, task))
    os.mkdir(parentfolder)

    # condorfile path
    condorfilename = os.path.join(parentfolder, "condor.submit")
    condorfile = open(condorfilename, "w")

    # Write Condor file
    print ("Executable   = " + exec_path, file = condorfile)
    print ("Universe     = Standard", file = condorfile)
    print ("InitialDir   = " + parentfolder, file = condorfile)
    # If a job is running for more than 10 minutes, kill it
    print ("periodic_remove = JobStatus == 2 && CurrentTime-EnteredCurrentStatus > 600", file = condorfile)
    print ("", file = condorfile)

    # Add Queue for each energy
    for i in range(0, n):
        # no verbose, relaxation
        arguments = "-d 0"
        # eos
        arguments += " -f " + eos_path
        # energy
        energy = i*step + e1
        arguments += " -e {:.2e}".format(energy)
        # Static sequence means that the ratio is 1
        arguments += " -r 1"
        print ("Arguments   = " + arguments, file = condorfile)
        print ("Output      = {}_{:.2e}.out".format(str(i).zfill(len(str(n))), energy),\
               file = condorfile)
        print ("Queue" , file = condorfile)
        print ("", file = condorfile)

    # Write a command.submit. condor.sumbit will be deleted to spare storage,
    # but it's important to know which command had been iussed to condor, this
    # is stored in command.submit
    commandfilename = os.path.join(parentfolder, "command.submit")
    commandfile = open(commandfilename, "w")
    print ("Arguments   = " + arguments, file = commandfile)

    condorfile.close()
    commandfile.close()

    if args.condor:
            # This is not the best wat to do that, but it works
            result = subprocess.check_output("condor_submit " + condorfilename, shell = True)
            # Delete condor.submit
            os.remove(condorfilename)
            # Print to STDIO results
            print(result.decode('utf8'))
    else:
        print ("condor.submit have been produced")

    # It is simpler to treat this case separatly and thus exit the program
    sys.exit(0)

# Set a to 0 if -R is not provided
if (not args.diff and not args.diff3):
    na  = 1
    na1 = 1
    na2 = 1
    nb  = 1
    a1initial = 0.
    a2initial = 0.
    a1final = 0.
    a2final = 0.
    b1 = 0.
    b2 = 0.

if (args.diff):
    nb  = 1
    na1 = na
    na2 = 1
    a1initial = ai
    a2initial = 0.
    a1final = af
    a2final = 0.
    b1 = 1.
    b2 = 1.

if (a1final == None): a1final = a1initial
if (a2final == None): a2final = a2initial
if (b2 == None): b2 = b1

# Compute B step
if (nb == 1):
    stepb = 0
else:
    stepb = (b2 - b1)/(nb - 1)

# Compute A1 step
if (na1 == 1):
    stepa1 = 0
else:
    stepa1 = (a1final - a1initial)/(na1 - 1)

# Compute A2 step
if (na2 == 1):
    stepa2 = 0
else:
    stepa2 = (a2final - a2initial)/(na2 - 1)

# Compute mass or jmoment step
if (N == 1):
    stepm = 0
    stepj = 0
else:
    if (task == "rmass"): stepm = (m2 - m1)/(N - 1)
    if (task == "jmoment"): stepj = (j2 - j1)/(N - 1)

# Compute energy step
if (n == 1):
    stepe = 0
else:
    stepe = (e2 - e1)/(n - 1)

# condorfiles will contain all the condor.submit paths produced
condorfiles = []

# basefolder summarizes how the script has been run
if (task == "rmass"):
    suffix = "_{}_Mi_{:.2}_Mf_{:.2}_ntot_{}/".format(task, m1, m2, na1*na2*nb)
if (task == "jmoment"):
    suffix = "_{}_Ji_{:.2}_Jf_{:.2}_ntot_{}/".format(task, j1, j2, na1*na2*nb)

fld = "{}_A1i_{:.2}_A1f_{:.2}_A2i_{:.2}_A2f_{:.2}_Bi_{:.2}_Bf_{:.2}_Ei_{:.2e}_Ef_{:.2e}".format(eos, a1initial, a1final, a2initial, a2final, b1, b2, e1, e2)  + suffix
basefolder = os.path.join(basefolder, fld)
os.mkdir(basefolder)

# Loop over values of A1
for q in range(0, na1):
    # Current value of A1
    a1 = q*stepa1 + a1initial

    fld = "{}_A1_{:.3}_A2i_{:.3}_A2f_{:.3}_Bi_{:.3}_Bf_{:.3}_Ei_{:.3e}_Ef_{:.3e}".format(eos, a1, a2initial, a2final, b1, b2, e1, e2) + suffix
    a1_folder = os.path.join(basefolder, fld)

    os.mkdir(a1_folder)
    # Loop over values of A2
    for l in range(0, na2):
        # Current value of A1
        a2 = l*stepa2 + a2initial
        fld = "{}_A1_{:.3}_A2_{:.3}_Bi_{:.3}_Bf_{:.3}_Ei_{:.3e}_Ef_{:.3e}".format(eos, a1, a2, b1, b2, e1, e2) + suffix
        a2_folder = os.path.join(basefolder, a1_folder, fld)
        os.mkdir(a2_folder)
        # Compute model with various value of beta
        for u in range(0,nb):
            # Current value of b
            b = u*stepb + b1
            fld = "{}_A1_{:.3}_A2_{:.3}_B_{:.3}_Ei_{:.3e}_Ef_{:.3e}".format(eos, a1, a2, b, e1, e2) + suffix
            b_folder = os.path.join(basefolder, a1_folder, a2_folder, fld)
            os.mkdir(b_folder)
            # loop over values of mass, or jmoment
            for k in range(0, N):
                if (task == "rmass"):
                    # Current value of mass
                    mass = k*stepm + m1
                    fld = "{}_A1_{:.3}_A2_{:.3}_B_{:.3}_Ei_{:.3e}_Ef_{:.3e}_M_{:.3}".format(eos, a1, a2, b, e1, e2, mass)
                    data_folder =  os.path.join(basefolder, a1_folder, a2_folder, b_folder, fld)

                if (task == "jmoment"):
                    # Current mass
                    jmoment = k*stepj + j1
                    fld = "{}_A1_{:.3}_A2_{:.3}_B_{:.3}_Ei_{:.3e}_Ef_{:.3e}_J_{:.3}".format(eos, a1, a2, b, e1, e2, jmoment)
                    data_folder =  os.path.join(basefolder, a1_folder, a2_folder, b_folder, fld)

                # Folder where are stored results
                os.mkdir(data_folder)

                # Produce condor.submit
                # condorfile path
                condorfilename = os.path.join(data_folder, "condor.submit")
                condorfiles.append(condorfilename)
                condorfile = open(condorfilename, "w")

                # Write Condor file
                print ("Executable   = " + exec_path, file = condorfile)
                print ("Universe     = Standard", file = condorfile)
                print ("InitialDir   = " + data_folder, file = condorfile)
                # If a job is running for more than 10 minutes, kill it
                print ("periodic_remove = JobStatus == 2 && CurrentTime-EnteredCurrentStatus > 600", file = condorfile)
                # print ("Notification = Complete", file = condorfile)
                # print ("notify_user  = spammozzola@gmail.com", file = condorfile)
                print ("", file = condorfile)

                for i in range(0, n):
                    # no verbose, relaxation
                    arguments = "-c 0.5 -d 0"
                    # eos
                    arguments += " -f " + eos_path
                    if (args.diff):  arguments += " -R diff -A {:.5}".format(a1)
                    if (args.diff3):  arguments += " -R diff -A {:.5} -D {:.5} -b {:.5}".format(a1, a2*a2, b)
                    # target
                    if (task == "rmass"):    arguments += " -s 1 -n 1 -p {:.5f}".format(mass)
                    if (task == "jmoment"):  arguments += " -s 2 -n 1 -p {:.5f}".format(jmoment)
                    # energy
                    energy = i*step + float(e1)
                    arguments += " -e {:.5e}".format(energy)
                    print ("Arguments   = " + arguments, file = condorfile)
                    print ("Output      = {}_{:.3e}.out".format(str(i).zfill(len(str(n))), energy),\
                           file = condorfile)
                    print ("Queue" , file = condorfile)
                    print ("", file = condorfile)

                # Write a command.submit. condor.sumbit will be deleted to spare storage,
                # but it's important to know which command had been iussed to condor, this
                # is stored in command.submit
                commandfilename = os.path.join(data_folder, "command.submit")
                commandfile = open(commandfilename, "w")
                print ("Arguments   = " + arguments, file = commandfile)

                condorfile.close()
                commandfile.close()

# It can take a while... 3 hours for 100k jobs
if args.condor:
    for f in condorfiles:
        # This is not the best wat to do that, but it works
        result = subprocess.check_output("condor_submit " + f, shell = True)
        # Delete condor.submit
        os.remove(f)
        # Print to STDIO results
        print(result.decode('utf8'))
        # Wait 10 s so that Conodor can handle better the queue
        time.sleep(10)
else:
    print ("condor.submit have been produced")
