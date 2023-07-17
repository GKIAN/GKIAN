#!/bin/env python3

import dispec as spec
import numpy as np
import matplotlib.pyplot as plt

spec.paramod.paragetarguments(input = '../data/input.nml',
  model = '../data/model.dat', deriv = True)
spec.math.mathinitialize()
spec.paramod.parainitialize()
f = spec.paramod.f.copy()
c = spec.paramod.c.copy()

spec.grtcmod.grtcinitialize()
spec.grtcmod.grtcevaluate()

def normlz(dat, deriv = False):
  if deriv:
    dat_ = np.log(np.abs(dat))
  else:
    dat_ = dat.copy()
    for i in range(dat.shape[0]):
      dat_[i, :] = dat[i, :] / np.amax(dat[i, :])
  return dat_

ker = normlz(spec.grtcmod.kermat)
plt.pcolormesh(f, c, ker.T, cmap = 'jet', shading = 'auto')
plt.xlabel('frequency (Hz)')
plt.ylabel('phase velocity (m/s)')
plt.colorbar()
plt.show()

if(spec.paramod.hasderiv):
  for i in range(4):
    dker = normlz(spec.grtcmod.dkermat[:, :, i], deriv = True)
    plt.pcolormesh(f, c, dker.T, cmap = 'jet', shading = 'auto')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('phase velocity (m/s)')
    plt.title(r'derivative at $ j = {:d} $'.format(i + 1))
    plt.colorbar()
    plt.show()

spec.grtcmod.grtcfinalize()
spec.paramod.parafinalize()
