#!/bin/bash

#/ Usage: simsync <sim_name> (eg RNS1)
#/ Copy form remote the output of the simulation <sim_name>
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

mkdir -p "$local_simulation_folder""/""$1"

watched_folder="$remote_simulation_folder""/$1"


is_ok=$(ssh "$user""@""$remote" "if [ -d '$watched_folder' ]; then echo '1'; fi;")

if [ -z "$is_ok" ]; then
   echo "Simulation does not exist on $(echo $remote | cut -d '.' -f 1)"
   exit 1
fi

if [ "$is_ok" -eq '1' ]; then
    folder=$(ssh "$user""@""$remote" "ls '$remote_simulation_folder''/''$1''/output-0000' | grep par | cut -d '.' -f 1")
    # h = human readable
    # v = verbose
    # P = partial (keep partially transferred files) and progress
    # z = compress
    # a = archive mode
    # e ssh = use ssh
    rsync -avzhPe ssh  "$user""@""$remote"":""$remote_simulation_folder""/""$1""/output-0000/""$folder""/data" \
          "$local_simulation_folder""/""$1"
fi
