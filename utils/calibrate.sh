#!/bin/bash
# Calibrate both camera and projector using checkerboards

if ! [ $# -eq 3 ]
then
    echo "usage: $(basename "$0") chk/ camChk.pgm projChk.pgm"
    echo
    echo -e "\toutputs: cam.dat proj.dat"
    echo -e '\tchk/ => checkerboards used to compute internal parameters'
    echo -e '\tchk/ should contain two subfolders : cam/ and proj/'

    exit 1
fi

utils_dir="$(dirname $(readlink -f $0))"

calibrate="$utils_dir/calibrate.py"

# Calibrate camera
$calibrate -e $2 $1/cam/* > cam.dat

# Calibrate projector
$calibrate -p -e $3 $1/proj/* > proj.dat
