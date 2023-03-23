#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

dat = np.loadtxt('gPS02.input.ker')

f = dat[1:, 0]
c = dat[0, 1:]
ker = dat[1:, 1:]

def normlz(dat):
  dat_ = dat.copy()
  for i in range(dat.shape[0]):
    dat_[i, :] = dat[i, :] / np.amax(dat[i, :])
  return dat_

ker = normlz(ker)

plt.pcolormesh(f, c, ker.T, cmap = 'jet', shading = 'auto')
plt.colorbar()
if len(sys.argv) > 1:
  plt.savefig(sys.argv[1] + '.png', dpi = 150)
else:
  plt.show()

# vim:ft=python tw=80 ts=4 sw=2 et ai
