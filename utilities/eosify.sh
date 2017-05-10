#!/bin/sh

#/ Usage: eosify <eos> (eg eosC)
#/ Produce a RNS 4.0 compatible eos from a RNS 1.1c EOS table
#/ Needs to be tweaked inside for paths
usage() {
    grep '^#/' "$0" | cut -c4-
    exit 0
}
expr "$*" : ".*--help" > /dev/null && usage

# Where the RNS 1.1c EOS are stored
eos_folder=/home/sbozzolo/master_thesis/rns4.0/eos_old/

# Where the RNS 4.0 table will be saved
out_folder=/home/sbozzolo/master_thesis/rns4.0/EOS/

# Path of HnG (executable found in RNS 4.0 code)
HnG=/home/sbozzolo/master_thesis/rns4.0/EOS/HnG/HnG

# Number of points
n=$(head -n 1 "$eos_folder""$1")

# Extract Pressure and Energy
awk '{print $1, $2}' "$folder""$1" | tail -n +2 > /tmp/tmp.$$

# Convert them into the new format
$HnG /tmp/tmp.$$ "$n" >> "$out_folder""$1"

# Remove tmp files
rm /tmp/tmp.$$
