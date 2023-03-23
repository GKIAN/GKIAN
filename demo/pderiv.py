#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

nlyr = 4
dat = np.loadtxt('gPS02.input.ker')
f = dat[1:, 0]
c = dat[0, 1:]

nf = len(f)
nc = len(c)

dg = np.fromfile('dgPS02.input.ker', dtype = np.double)
dg = dg.reshape((nlyr, nc, nf))

def normlz(dat):
  dat_ = np.log(np.abs(dat))
  return dat_

for i in range(nlyr):
  dp = normlz(dg[i, :, :])
  plt.pcolormesh(f, c, dp, cmap = 'jet', shading = 'auto')
  plt.colorbar()
  plt.xlabel('frequency (Hz)')
  plt.ylabel('phase velocity (m/s)')
  plt.title(r'derivative at $ j = {:d} $'.format(i + 1))
  if len(sys.argv) > 1:
    plt.savefig('{:s}-dgPS02-{:d}.png'.format(sys.argv[1], i + 1), dpi = 150)
    plt.clf()
  else:
    plt.show()

# vim:ft=python tw=80 ts=4 sw=2 et ai
