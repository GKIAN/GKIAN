#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Mon 12 Dec 2022 09:22:49 PM CST
#-------------------------------------------------------------------------------

import dispec as spec
import numpy as np
import matplotlib.pyplot as plt
import argparse, json
import os, sys

dparfile = 'defpars/senker.json'
if not os.path.isfile(dparfile):
  pyhome = os.path.split(sys.argv[0])[0]
  dparfile = pyhome + '/../defpars/senker.json'
with open(dparfile) as fin:
  dpars = json.load(fin)

#===============================================================================

ap = argparse.ArgumentParser(description = 'to calculate and plot spectrum ' +
  'sensitivity of an input model.', prefix_chars = '-+')
# for 'readfile' group
ag1 = ap.add_argument_group('readfile', 'read from a sen-ker file to plot.')
ag1.add_argument('+k', '--kerfile', metavar = 'senKernelFile',
  default = dpars['+k'], help = 'read from the sen-ker file; otherwise, ' +
  're-calculate spectrum sensitivity for an input model. ' +
  '[default: %(default)s]')
# for 'evaluate' group
ag2 = ap.add_argument_group('evaluate', 'evaluate for the input model to plot.')
ag2.add_argument('modfile', metavar = 'inputModelFile', nargs = '?',
  default = dpars['EP'], help = 'the input model file. [default: %(default)s]')
aeg = ag2.add_mutually_exclusive_group()
aeg.add_argument('-I', '--inml', metavar = 'inputNMLFile',
  default = dpars['-I'], help = 'the input NAMELIST file. ' +
  '[default: %(default)s]')
aeg.add_argument('-W', '--ifactor', metavar = 'imgFreqFactor',
  default = dpars['-W'], type = float, help = 'the imaginary frequency ' +
  'factor. If not None, read f/c sampling points from the respective ' +
  'input files. [default: %(default)s]')
ag2.add_argument('-f', '--ffile', metavar = 'freqFile', default = dpars['-f'],
  help = 'the frequency input file. [default: %(default)s]')
ag2.add_argument('-c', '--cfile', metavar = 'cFile', default = dpars['-c'],
  help = 'the phase velocity input file. [default: %(default)s]')
ag2.add_argument('-i', '--fidxargs', metavar = 'iStart[:iEnd[:iStride]]',
  default = dpars['-i'], help = 'the arguments for calculating the index of ' +
  'frequency sampling-points in inversion. [default: %(default)s]')
ag2.add_argument('-j', '--cidxargs', metavar = 'iStart[:iEnd[:iStride]]',
  default = dpars['-j'], help = 'the arguments for calculating the index of ' +
  'phase velocity sampling-points in inversion. [default: %(default)s]')
ag2.add_argument('-w', '--winargs', default = dpars['-w'],
  metavar = 'iStart[:iEnd[:iStride]],jStart[:jEnd[:jStride]]',
  help = 'the arguments for windowing spectrogram in inversion. ' +
  '[default: %(default)s]')
ag2.add_argument('-u', '--notsiu', action = 'store_true',
  help = 'NOT SI Units used in the initial model, but common units, such as ' +
  'km, km/s, g/cm^3.')
ag2.add_argument('+K', '--Kerfile', metavar = 'senKernelFile',
  default = dpars['+K'], help = 'save evaluated spectrum sensitivity result ' +
  'to the file; otherwise, not save. [default: %(default)s]')
# for extracting
ap.add_argument('-E', '--extkey', metavar = 'key=Value', default = dpars['-E'],
  help = 'the profile parameters for spectrum sensitivity. Acceptable key: ' +
  'f, c and z. [default: %(default)s]')
# for plotting
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-s', '--savepath', metavar = 'saveFigPath',
  default = dpars['-s'], help = 'save figure(s) to the path. ' +
  '[default: %(default)s]')
ag.add_argument('-S', '--show', action = 'store_true',
  help = 'show figure(s) immediately.')
ap.add_argument('-O', '--outpref', metavar = 'outFigPrefix',
  default = dpars['-O'], help = 'save figure as a file with ' +
  'the name prefix. [default: %(default)s]')
ap.add_argument('-p', '--dpi', metavar = 'dpiValue', type = int,
  default = dpars['-p'], help = 'specify the value of dpi for saving figure. ' +
  '[default: %(default)s]')
ap.add_argument('-k', '--kmunit', action = 'store_true',
  help = 'use km as the length unit.')
ap.add_argument('-e', '--epsfmt', action = 'store_true',
  help = 'save figure(s) with EPS format; otherwise, with PNG format.')
args = ap.parse_args()

if args.kmunit:
  lu2p = 1.0e-3
  lunit = 'km'
else:
  lu2p = 1.0
  lunit = 'm'

if args.epsfmt:
  psfx = '.eps'
else:
  psfx = '.png'

#===============================================================================

def colon2lint(string):
  subs = string.split(':')
  ret = [ 0, 0, 1 ]
  for i in range(len(subs)):
    ret[i] = int(subs[i])
  return ret

if args.kerfile is None:
  spec.math.mathinitialize()
  if args.ifactor is None:
    spec.paramod.paragetarguments(input = args.inml, model = args.modfile,
      deriv = True)
    spec.paramod.parainitialize()
  else:
    spec.paramod.paragetarguments(model = args.modfile, deriv = True)
    spec.paramod.parainitialize(args.ifactor)
    spec.paramod.f = np.loadtxt(args.ffile)
    spec.paramod.fnum = len(spec.paramod.f)
    spec.paramod.c = np.loadtxt(args.cfile)
    spec.paramod.cnum = len(spec.paramod.c)

  if args.notsiu:
    u2si = 1.0e+3
  else:
    u2si = 1.0
  spec.paramod.z *= u2si
  spec.paramod.h *= u2si
  spec.paramod.c *= u2si

  if args.winargs is None:
    fias = colon2lint(args.fidxargs)
    cias = colon2lint(args.cidxargs)
  else:
    winargs = args.winargs.split(',')
    fias = colon2lint(winargs[0])
    cias = colon2lint(winargs[1])

  if fias[1] <= 0:
    fias[1] += spec.paramod.fnum
  if cias[1] <= 0:
    cias[1] += spec.paramod.cnum

  fs = spec.paramod.f.copy()
  cs = spec.paramod.c.copy()
  z  = spec.paramod.z.copy()
  z[-1] = z[-2] * 2.0 - z[-3]

  f = fs[ fias[0]:fias[1]:fias[2] ]
  c = cs[ cias[0]:cias[1]:cias[2] ]
  spec.paramod.f = f.copy()
  spec.paramod.fnum = len(f)
  spec.paramod.c = c.copy()
  spec.paramod.cnum = len(c)

  spec.grtcmod.grtcinitialize()
  spec.grtcmod.grtcevaluate()
  dker = spec.grtcmod.dkermat.copy()

  if not args.Kerfile is None:
    np.savez(args.Kerfile, u = u2si, f = f, c = c, z = z, ker = dker)

  spec.grtcmod.grtcfinalize()
  spec.paramod.parafinalize()

else:
  data = np.load(args.kerfile)
  u2si = data['u']
  f = data['f']
  c = data['c']
  z = data['z']
  dker = data['ker']

#===============================================================================

labels = { 'f': 'Frequency (Hz)', 'c': 'Phase velocity (%s/s)' % (lunit),
  'z': 'Depth (%s)' % (lunit) }
yinv = True

if args.extkey is None:
  vker = ( np.sum( np.abs(dker), axis = 1 ) ).squeeze()
  #vker = ( np.sum( dker, axis = 1 ) ).squeeze()
  x = f
  y = z * lu2p
  xlabel = labels['f']
  ylabel = labels['z']
  title = 'The sum of absolute value among phase-velocity axis'
else:
  kvstr = args.extkey.split('=')
  kext = kvstr[0]
  vext = float(kvstr[1])
  if kext == 'f':
    fi = np.argmin( np.abs(f - vext) )
    vker = dker[fi, :, :]
    x = c * lu2p
    y = z * lu2p
    xlabel = labels['c']
    ylabel = labels['z']
    title = 'The profile for f = %.3f Hz' % (f[fi])
  else:
    vext *= u2si
    if kext == 'c':
      ci = np.argmin( np.abs(c - vext) )
      vker = dker[:, ci, :]
      x = f
      y = z * lu2p
      xlabel = labels['f']
      ylabel = labels['z']
      title = 'The profile for c = %.3f %s/s' % (c[ci] * lu2p, lunit)
    elif kext == 'z':
      yinv = False
      if vext > z[-2]:
        zi = -2
      else:
        zi = np.searchsorted(z, vext, side = 'right') - 1
      vker = dker[:, :, zi]
      x = f
      y = c * lu2p
      xlabel = labels['f']
      ylabel = labels['c']
      title = 'The profile for z in [%.3f, %.3f] %s' % (z[zi] * lu2p,
        z[zi + 1] * lu2p, lunit)

#===============================================================================

vmax = np.amax(np.abs(vker))
plt.pcolormesh(x, y, vker.T, cmap = 'seismic', shading = 'auto',
  rasterized = True)
plt.colorbar()
plt.clim(( - vmax, vmax))
if yinv:
  plt.gca().invert_yaxis()
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
if args.show:
  plt.show()
else:
  plt.savefig(args.savepath + args.outpref + psfx, dpi = args.dpi,
    bbox_inches = 'tight')

# vim:ft=python tw=80 ts=4 sw=2 et ai
