#!/usr/bin/env python3

# Condor Parser for RNS v4 and HTCondor 8.4.5
# Prepare a bunch of condor.sumbit and than submit them

# Author: Gabriele Bozzola (sbozzolo)
# Email: sbozzolator@gmail.com
# Version: 1.5
# First Stable: 13/03/17
# Last Edit: 23/03/17

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
parser.add_argument("-t", "--task", type = str, choices = ["static", "kepler", "rmass", "jmoment"],
                    required = True, help = "define which kind of sequences have to be computed")
parser.add_argument("-o", "--out", type = str,
                    help = "set output folder")
parser.add_argument('-N', "--nsequences", action="store", type = int,
                    help =
                    "set number of sequences")
parser.add_argument('-n', "--nmodels", action="store", type = int,
                    required = True, help = "set number of models for each sequence")
parser.add_argument('-na', "--namodels", action="store", type = int,
                    help = "set number of models for each sequence of rotation")
parser.add_argument('-e1', "--energy1", action="store", type = float,
                    required = True, help = "initial central energy")
parser.add_argument('-e2', "--energy2", action="store", type = float,
                    required = True, help = "final central energy")
parser.add_argument("-R", "--diff", help = "enable differential rotation",
                    action = "store_true")
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
parser.add_argument("-c", "--condor", help = "sumbit jobs to Condor",
                    action = "store_true")
parser.add_argument('--version', action='version', version='%(prog)s 1.5')

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
a1   = args.rotation1
a2   = args.rotation2

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

# Remove old results, if they exist
if os.path.exists(basefolder):
    shutil.rmtree(basefolder)
os.mkdir(basefolder)

# Set a to 0 if -R is not provided
if (not args.diff):
    na = 1
    stepa = 0
    a1 = 0.

# Compute A step
if (na == 1):
    stepa = 0
else:
    stepa = (a2 - a1)/(na - 1)

# condorfiles will contain all the condor.submit paths produced
condorfiles = []

# Loop over values of A
for j in range(0, na):

    # Current value of a
    a = j*stepa + a1

    # There's no need to compute more than 1 models if the task id static or kepler
    if (task == "static" or task == "kepler"): N = 1

    # parentfolder is the folder with fixed A
    # parantfodlder is customized depending on the task, so it can be read in other scripts
    if (task == "rmass"):
        parentfolder = os.path.join(basefolder,
                                  "{}_n{}_A_{:.2}_Ei{:.1}_Ef{:.1}_{}_Mi{}_Mf{}/".format(eos, n, a, e1, e2, task, m1, m2))
    if (task == "jmoment"):
        parentfolder = os.path.join(basefolder, \
                                    "{}_n{}_A_{:.2}_Ei{:.1}_Ef{:.1}_{}_Ji{}_Jf{}/".format(eos, n, a, e1, e2, task, j1, j2))
    if (task == "static"):
        parentfolder = os.path.join(basefolder, \
                                    "{}_n{}_Ei{:.1}_Ef{:.1}_{}/".format(eos, n, e1, e2, task))
    if (task == "kepler"):
        if (args.diff):
            parentfolder = os.path.join(basefolder, \
                                       "{}_n{}_A_{:.2}_Ei{:.1}_Ef{:.1}_{}/".format(eos, n, a, e1, e2, task))
        else:
            parentfolder = os.path.join(basefolder, \
                                    "{}_n{}_Ei{:.1}_Ef{:.1}_{}/".format(eos, n, e1, e2, task))

    os.mkdir(parentfolder)

    # loop over values of mass, or jmoment
    for k in range(0, N):

        if (task == "rmass"):
            # Mass steps
            if (N == 1):
                # if N = 1 there's no need to step
                stepm = 0
            else:
                stepm = (m2 - m1)/(N - 1)

        if (task == "jmoment"):
            # Moment steps
            if (N == 1):
                # if N = 1 there's no need to step
                stepj = 0
            else:
                stepj = (j2 - j1)/(N - 1)

        # Energy step
        step = (e2 - e1)/(n - 1)

        # folder is the name of the folder which will contain actual data

        if (task == "rmass"):
            # Current mass
            mass = k*stepm + m1
            folder = "{}_n{}_A_{:.2}_Ei{:.1}_Ef{:.1}_{}_M{:.2f}".format(eos, n, a, e1, e2, task, mass)

        if (task == "jmoment"):
            # Current mass
            jmoment = k*stepj + j1
            folder = "{}_n{}_A_{:.2}_Ei{:.1}_Ef{:.1}_{}_J{:.2f}".format(eos, n, a, e1, e2, task, jmoment)

        if (task == "static"):
            folder = "{}_n{}_A_{:.2}_Ei{:.1}_Ef{:.1}_{}".format(eos, n, a, e1, e2, task)

        if (task == "kepler"):
            folder = "{}_n{}_Ei{:.1}_Ef{:.1}_{}".format(eos, n, a, e1, e2, task)
            if (args.diff): folder = "/{}_n{}_A_{:.2}_Ei{:.1}_Ef{:.1}_{}/".format(eos, n, a, e1, e2, task)

        os.mkdir(os.path.join(parentfolder, folder))

        # condorfile path
        condorfilename = os.path.join(parentfolder, folder, "condor.submit")
        condorfiles.append(condorfilename)
        condorfile = open(condorfilename, "w")

        # Write Condor file
        print ("Executable   = " + exec_path, file = condorfile)
        print ("Universe     = Standard", file = condorfile)
        print ("InitialDir   = " + os.path.join(parentfolder, folder), file = condorfile)
        # If a job is running for more than 10 minutes, kill it
        print ("periodic_remove = JobStatus == 2 && CurrentTime-EnteredCurrentStatus > 300", file = condorfile)
        # print ("Notification = Complete", file = condorfile)
        # print ("notify_user  = spammozzola@gmail.com", file = condorfile)
        print ("", file = condorfile)

        for i in range(0, n):
            # no verbose, relaxation
            arguments = "-c 0.5 -d 0"
            # eos
            arguments += " -f " + eos_path
            if (args.diff):  arguments += " -R diff -A {:.2}".format(a)
            # target
            if (task == "rmass"):    arguments += " -s 1 -n 1 -p {:.2f}".format(mass)
            if (task == "jmoment"):  arguments += " -s 2 -n 1 -p {:.2f}".format(jmoment)
            if (task == "kepler"):   arguments += " -s 3 -n 1"
            if (task == "static"):   arguments += " -r 1"
            # energy
            energy = i*step + float(e1)
            arguments += " -e {:.2e}".format(energy)
            print ("Arguments   = " + arguments, file = condorfile)
            print ("Output      = {}_{:.2e}.out".format(str(i).zfill(len(str(n))), energy),\
                   file = condorfile)
            print ("Queue" , file = condorfile)
            print ("", file = condorfile)

        # Write a command.submit. condor.sumbit will be deleted to spare storage,
        # but it's important to know which command had been iussed to condor, this
        # is stored in command.submit
        commandfilename = os.path.join(parentfolder, folder, "command.submit")
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
