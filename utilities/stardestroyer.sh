#!/bin/bash

#/ Usage: stardestroyer
#/ Delete every occurence of the folder Output_star to save storage
#/ Should always run in background
usage() {
    grep '^#/' "$0" | cut -c4-
    exit 0
}
expr "$*" : ".*--help" > /dev/null && usage

while :
do
    find . -name "Output_star" -type d -exec rm -rf {} +
    echo "Destroyed! `date`"
    sleep 900
done
