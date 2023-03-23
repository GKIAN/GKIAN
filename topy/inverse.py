
#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Tue 27 Apr 2021 04:42:32 PM CST
#-------------------------------------------------------------------------------

import empirical as emp
import numpy as np
import scipy.linalg as sla
import os

class Inversion:
  def __init__(self, invargs):
    (spec, cmdargs, nw, fa, ca, wwin, data, xwlev, whoami) = invargs
    self.itr = 0
    self.spec = spec
    self.nlyr = spec.paramod.nlayer
    self.fs = spec.paramod.f.copy()
    self.cs = spec.paramod.c.copy()
    self.d2init = cmdargs.regargs[0]
    self.smooth = cmdargs.regargs[1]
    self.outpath = cmdargs.outpath
    self.verbosity = cmdargs.verbosity
    self.nottn = cmdargs.nottn
    self.emptype = cmdargs.emptype
    self.afwgt = cmdargs.fwgtarg
    self.nw = nw
    self.fa = fa
    self.ca = ca
    self.wwin = wwin
    self.data = data
    self.xwlev = xwlev
    self.whoami = whoami
    self.model = np.zeros((self.nlyr, 5))
    self.model[:, 0] = np.arange(self.nlyr) + 1
    self.model[:, 1] = spec.paramod.z[: - 1]
    self.model[:, 2] = spec.paramod.rho.copy()
    self.model[:, 3] = spec.paramod.beta.copy()
    self.model[:, 4] = spec.paramod.alpha.copy()
    self.xbase = spec.paramod.beta.copy()
    self.xmisf = np.nan * np.ones(self.nlyr)
    self.gker = [ None ] * self.nw
    self.dker = [ None ] * self.nw
    self.maxk = [ None ] * self.nw
    self.nfc = [ self.data[iw].size for iw in range(self.nw) ]
    tv = np.zeros(self.nlyr)
    tv[0] = 0.5; tv[1] = - 0.25
    self.tplz = sla.toeplitz(tv)[1: - 1, :]
    self.isroot = (self.whoami == 0)
    for iw in range(self.nw):
      self.normalize(self.data[iw])

  def prepare(self, index, xinit = None):
    self.itr = 0
    if not xinit is None:
      self.update_model(xinit)
    if self.verbosity > 0:
      if self.isroot:
        print('Now running for the %d-th model ...' % (index), flush = True)
    if self.verbosity > 1:
      self.subpath = self.outpath + 'i%03d/' % (index)
      os.makedirs(self.subpath, exist_ok = True)
      self.save_model(self.subpath + 'initmodel.txt')

  def calc_misfit(self, x):
    return np.nan

  def calc_gradient(self, x):
    return np.nan * np.ones(self.nlyr)

  def update_model(self, x, isget = False):
    beta = self.xbase + self.xwlev * x
    if isget:
      retv = beta.copy()
      return retv
    else:
      if self.emptype == 'B':
        self.model[:, 2:5] = emp.brocherFit(beta)
      elif self.emptype == 'U':
        self.model[:, 2:5] = emp.userFit(beta)
      else:
        self.model[:, 3] = beta
      self.spec.paramod.rho   = self.model[:, 2]
      self.spec.paramod.beta  = self.model[:, 3]
      self.spec.paramod.alpha = self.model[:, 4]
      return None

  def save_model(self, outputfile):
    fout = open(outputfile, 'wt')
    fout.write('Layer-No.   Depth(m)   Densities(kg/m3)  S-Velocity(m/s)  ' +
      'P-Velocity(m/s)\n')
    np.savetxt(fout, self.model, fmt = '%5d      %10.2f   %14.3f    ' +
      '%13.3f    %13.3f')
    fout.close()

  def update_window(self, iw):
    self.spec.paramod.fnum = self.data[iw].shape[0]
    self.spec.paramod.cnum = self.data[iw].shape[1]
    self.spec.paramod.f = self.fs[ self.fa[iw][0]:self.fa[iw][1]:self.fa[iw][2] ]
    self.spec.paramod.c = self.cs[ self.ca[iw][0]:self.ca[iw][1]:self.ca[iw][2] ]
    self.spec.grtcmod.kermat = np.zeros_like(self.data[iw])
    self.spec.grtcmod.dkermat = np.zeros( (self.data[iw].shape[0],
      self.data[iw].shape[1], self.nlyr) )

  def normalize(self, what):
    fctr = np.ones_like(what)
    if self.nottn:
      for i in range(what.shape[0]):
        wmin = np.amin(what[i, :])
        wmax = np.amax(what[i, :])
        fctr[i, :] = wmax - wmin
        what[i, :] = (what[i, :] - wmin) / fctr[i, :]
    else:
      what[:, :] = np.where(what < 0.0, 0.0, what)
      for i in range(what.shape[0]):
        fctr[i, :] = np.amax(np.abs(what[i, :]))
        what[i, :] /= fctr[i, :]
    return fctr

  def regularizate(self, x, isderiv = False, pure = 'A'):
    # b = self.update_model(x, True)
    b = self.xwlev * x #TODO
    if not isderiv:
      # initiating regularizate
      regi = 0.5 * np.sum(x * x) / self.nlyr
      # smoothness regularization
      regs = 0.5 * np.sum((self.tplz @ b) ** 2) / (self.nlyr - 2)
      if pure == 'I':
        return regi
      elif pure == 'S':
        return regs
      else:
        return self.d2init * regi + self.smooth * regs
    else:
      regi = self.d2init * x / self.nlyr
      regs = self.smooth * ( (self.tplz @ b).T @ self.tplz ) \
        / (self.nlyr - 2)
      return regi + regs

  def callback(self, x):
    self.itr = self.itr + 1
    self.update_model(x)
    # print result to stdout
    if (self.verbosity > 0) & self.isroot:
      print('x = ', x)
      print('>> in iteration %4d: beta = ' % (self.itr),
        self.spec.paramod.beta, flush = True)
    # write result to file
    if self.verbosity > 1:
      self.save_model(self.subpath + 'invmodel-%03d.dat' % (self.itr))

#===============================================================================

class InversionSD2N(Inversion):
  # SD2N: Spectrum Difference's 2-Norm
  def __init__(self, invargs):
    super().__init__(invargs)
    self.diff = [ None ] * self.nw
    if self.afwgt is None:
      self.fwgt = [ np.ones_like(self.data[iw]) for iw in range(self.nw) ]
    else:
      # different weights for different frequencies
      self.fwgt = [ np.zeros_like(self.data[iw]) for iw in range(self.nw) ]
      for iw in range(self.nw):
        f = self.fs[ self.fa[iw][0]:self.fa[iw][1]:self.fa[iw][2] ]
        for i in range(self.data[iw].shape[0]):
          self.fwgt[iw][i, :] = ( 1.0 / f[i] ) ** self.afwgt

  def calc_misfit(self, x, usereg = True):
    self.xmisf = x.copy()
    self.update_model(x)
    self.misfit = 0.0
    for iw in range(self.nw):
      self.update_window(iw)
      self.spec.grtcmod.grtcevaluate()
      self.gker[iw] = self.spec.grtcmod.kermat.copy()
      self.dker[iw] = self.spec.grtcmod.dkermat.copy()
      self.maxk[iw] = self.normalize(self.gker[iw])
      self.diff[iw] = self.gker[iw] - self.data[iw]
      self.misfit += 0.5 * np.sum(self.diff[iw] * self.diff[iw] \
        * self.fwgt[iw]) / self.nfc[iw] * self.wwin[iw]
    if usereg:
      self.misfit += self.regularizate(x)
    return self.misfit

  def calc_gradient(self, x):
    self.grad = np.zeros(self.nlyr)
    if not np.allclose(x, self.xmisf):
      self.calc_misfit(x)
    for iw in range(self.nw):
      for i in range(self.nlyr):
        if self.nottn:
          dkeri = self.dker[iw][:, :, i]
        else:
          dkeri = np.where(self.gker[iw] < 0.0, 0.0, self.dker[iw][:, :, i]) #TODO
        ddk = self.diff[iw] * dkeri / self.maxk[iw]
        self.grad[i] += np.sum(ddk * self.fwgt[iw]) / self.nfc[iw] \
          * self.wwin[iw]
    self.grad += self.regularizate(x, isderiv = True)
    self.grad *= self.xwlev
    return self.grad

#-------------------------------------------------------------------------------

class InversionGCCC(Inversion):
  # GCCC: Global Cross-Correlation Coefficient
  def __init__(self, invargs):
    super().__init__(invargs)
    if not self.afwgt is None:
      for iw in range(self.nw):
        f = self.fs[ self.fa[iw][0]:self.fa[iw][1]:self.fa[iw][2] ]
        for i in range(self.data[iw].shape[0]):
          self.data[iw][i, :] = self.data[iw][i, :] ** (self.afwgt / f[i])
    self.ccd = [ np.sqrt(np.sum(self.data[iw] * self.data[iw]))
      for iw in range(self.nw) ]
    self.cck = np.zeros(self.nw)
    self.wcc = np.zeros(self.nw)

  def calc_misfit(self, x, usereg = True):
    self.xmisf = x.copy()
    self.update_model(x)
    self.misfit = 0.0
    for iw in range(self.nw):
      self.update_window(iw)
      self.spec.grtcmod.grtcevaluate()
      self.gker[iw] = self.spec.grtcmod.kermat.copy()
      self.dker[iw] = self.spec.grtcmod.dkermat.copy()
      self.maxk[iw] = self.normalize(self.gker[iw])
      self.cck[iw] = np.sqrt(np.sum(self.gker[iw] * self.gker[iw]))
      self.wcc[iw] = np.sum(self.gker[iw] * self.data[iw]) \
        / (self.cck[iw] * self.ccd[iw])
      self.misfit += (1.0 - self.wcc[iw]) * self.wwin[iw]
    if usereg:
      self.misfit += self.regularizate(x)
    return self.misfit

  def calc_gradient(self, x):
    self.grad = np.zeros(self.nlyr)
    if not np.allclose(x, self.xmisf):
      self.calc_misfit(x)
    for iw in range(self.nw):
      for i in range(self.nlyr):
        kdk = self.gker[iw] * self.dker[iw][:, :, i] / self.maxk[iw]
        ddk = self.data[iw] * self.dker[iw][:, :, i] / self.maxk[iw]
        self.grad[i] += ( np.sum(kdk) / (self.cck[iw] ** 2) * self.wcc[iw] \
          - np.sum(ddk) / (self.cck[iw] * self.ccd[iw]) ) \
          * self.wwin[iw]
    self.grad += self.regularizate(x, isderiv = True)
    self.grad *= self.xwlev
    return self.grad

# vim:ft=python tw=80 ts=4 sw=2 et ai
