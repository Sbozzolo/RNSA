#!/bin/sh

#/ Usage: rhoextractor <simulation>
#/ Produce a <simulation>.out file that contains the evolution of rho(t)/rho(0)
#/ Needs to be tweaked inside for paths
usage() {
    grep '^#/' "$0" | cut -c4-
    exit 0
}
expr "$*" : ".*--help" > /dev/null && usage

if [ "$#" -ne 1 ]; then
    echo "Illegal number of arugments. Try --help"
    exit 1
fi

# Where the are simulations?
sims=/home/sbozzolo/master_thesis/simulations

if [ -d "$sim""/""$i" ]; then
   echo "Simulation does not exist!"
   exit 1
fi

find "$sim""/""$i" -type f -name "hydrobase-rho.maximum.asc" -exec echo {} \; | sort | xargs cat | sed 's/#.*$//' | sed '/^$/d' | awk '{if (NR == 1) mass = $3; print $2/204, $3/mass}' > $i.dat
