#!/bin/bash

set -e
set -x

export OMP_NUM_THREADS=4

# ../../bin/grtcKer ../data/input.nml ../data/model.dat
../../bin/grtcKer -d ../data/input.nml ../data/model.dat

