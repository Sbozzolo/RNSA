#!/bin/bash

#/ Usage: watchsim <sim_name> (eg RNS1)
#/ Prints the log file
#/ remote_config needs to be tweaked inside for paths
usage() {
    grep '^#/' "$0" | cut -c4-
    exit 0
}
expr "$*" : ".*--help" > /dev/null && usage

source remote_config

if [ "$#" -ne 1 ]; then
    echo "Illegal number of arugments. Try --help"
    exit 1
fi

watched_folder="$remote_simulation_folder""/$1"

is_ok=$(ssh "$user""@""$remote" "if [ -d '$watched_folder' ]; then echo '1'; fi;")

if [ -z "$is_ok" ]; then
   echo "Simulation does not exist on $(echo $remote | cut -d '.' -f 1)"
   exit 1
fi

if [ "$is_ok" -eq '1' ]; then
    log_file="$watched_folder""/$1"".log"
    ssh "$user""@""$remote" "tail -f $log_file"
fi
