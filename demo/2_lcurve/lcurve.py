#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Sun 11 Dec 2022 04:55:06 PM CST
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as pl

dat = np.loadtxt('output/allinone.txt', skiprows = 1)[:3, :]

a = dat[0, :]
res = dat[1, :]
reg = dat[2, :]
na = len(a)

if 0:
  res = np.log(res)
  #reg = np.log(reg)
  reg = np.log(reg / 1000.0)

pl.plot(res, reg, 'o-')
for i in range(na):
  pl.text(res[i], reg[i], '%g' % (a[i]))
pl.show()

# vim:ft=python tw=80 ts=4 sw=2 et ai
