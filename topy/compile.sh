#!/bin/bash
f2py -m dispec -h temp.pyf src/constant.F90 src/math.F90 src/parameter.F90 src/grtc.F90
f2py -m dispec -c temp.pyf src/constant.F90 src/math.F90 src/parameter.F90 src/grtc.F90 --f90flags='-fopenmp' -lgomp -DIEEE -DDERIV=_BETA -DCOLORPRINT
