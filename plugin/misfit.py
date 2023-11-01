#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Fri 16 Apr 2021 16:14:52 PM CST
#===============================================================================

import dispec as spec
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import argparse, json
import glob, os, sys

dparfile = 'defpars/misfit.json'
if not os.path.isfile(dparfile):
  pyhome = os.path.split(sys.argv[0])[0]
  dparfile = pyhome + '/../defpars/misfit.json'
with open(dparfile) as fin:
  dpars = json.load(fin)

#===============================================================================

ap = argparse.ArgumentParser(description = 'to plot misfit.')
# for 'readfile' group
ag1 = ap.add_argument_group('readfile', 'read from a misfit file to plot.')
ag1.add_argument('-m', '--misfile', metavar = 'misfitFile',
  default = dpars['-m'], help = 'read from the misfit file; otherwise, ' +
  're-calculate misfit for all models. [default: %(default)s]')
# for 'evaluate' group
ag2 = ap.add_argument_group('evaluate', 'evaluate from the inversed results ' +
  'to plot.')
ag2.add_argument('fdata', nargs = '?', metavar = 'dataFile',
  default = dpars['EP'], help = 'the data file. [default: %(default)s]')
ag2.add_argument('-0', '--initial', metavar = 'initModel',
  default = dpars['-0'], help = 'the initial model. [default: %(default)s]')
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
ag2.add_argument('-R', '--outdir', metavar = 'outputDir', default = dpars['-R'],
  help = 'to plot all models under the inversion-output directory, ' +
  'recursively. [default: %(default)s]')
ag2.add_argument('-i', '--fidxargs', metavar = 'iStart[:iEnd[:iStride]]',
  default = dpars['-i'], help = 'the arguments for calculating the index of ' +
  'frequency sampling-points in inversion. [default: %(default)s]')
ag2.add_argument('-j', '--cidxargs', metavar = 'iStart[:iEnd[:iStride]]',
  default = dpars['-j'], help = 'the arguments for calculating the index of ' +
  'phase velocity sampling-points in inversion. [default: %(default)s]')
ag2.add_argument('-w', '--winargs', nargs = '+',
  metavar = 'iStart[:iEnd[:iStride]],jStart[:jEnd[:jStride]][,weight]',
  default = dpars['-w'], help = 'the arguments for windowing spectrogram ' +
  'in inversion. [default: %(default)s]')
ag2.add_argument('-t', '--datatp', action = 'store_true',
  help = 'transpose the data matrix.')
ag2.add_argument('-u', '--notsiu', action = 'store_true',
  help = 'NOT SI Units used in the initial model, but common units, such as ' +
  'km, km/s, g/cm^3.')
ag2.add_argument('-M', '--Misfile', metavar = 'misfitFile',
  default = dpars['-M'], help = 'save evaluated misfit result to the file; ' +
  'otherwise, not save. [default: %(default)s]')
# for plotting
ap.add_argument('-l', '--logy', action = 'store_true',
  help = "plot with 'semilogy'; otherwise, with 'plot'.")
aeg = ap.add_mutually_exclusive_group()
aeg.add_argument('-s', '--savepath', metavar = 'saveFigPath',
  default = dpars['-s'], help = 'save figure(s) to the path. ' +
  '[default: %(default)s]')
aeg.add_argument('-S', '--show', action = 'store_true',
  help = 'show figure(s) immediately.')
ap.add_argument('-O', '--outpref', metavar = 'outFigPrefix',
  default = dpars['-O'], help = 'save figure as a file with ' +
  'the name prefix. [default: %(default)s]')
ap.add_argument('-p', '--dpi', metavar = 'dpiValue', type = int,
  default = dpars['-p'], help = 'specify the value of dpi for saving figure. ' +
  '[default: %(default)s]')
ap.add_argument('-e', '--epsfmt', action = 'store_true',
  help = 'save figure(s) with EPS format; otherwise, with PNG format.')
args = ap.parse_args()

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

def tonorm(dk, nottn = False):
  if nottn:
    ndk = dk.copy()
    for i in range(ndk.shape[0]):
      dkmin = np.amin(ndk[i, :])
      dkmax = np.amax(ndk[i, :])
      ndk[i, :] = (ndk[i, :] - dkmin) / (dkmax - dkmin)
  else:
    ndk = np.where(dk < 0.0, 0.0, dk)
    for i in range(ndk.shape[0]):
      ndk[i, :] = ndk[i, :] / np.amax(np.abs(ndk[i, :]))
  return ndk

if args.misfile is None:
  if args.fdata[-2:] == 'h5':
    h5id = h5py.File(args.fdata, 'r')
    data = h5id['amp'][:]
    h5id.close()
  else:
    data = np.loadtxt(args.fdata)
  if args.datatp:
    data = data.T

  spec.math.mathinitialize()
  if args.ifactor is None:
    spec.paramod.paragetarguments(input = args.inml, model = args.initial)
    spec.paramod.parainitialize()
  else:
    spec.paramod.paragetarguments(model = args.initial)
    spec.paramod.parainitialize(args.ifactor)
    spec.paramod.f = np.loadtxt(args.ffile)
    spec.paramod.fnum = len(spec.paramod.f)
    spec.paramod.c = np.loadtxt(args.cfile)
    spec.paramod.cnum = len(spec.paramod.c)

  if args.notsiu:
    u2si = 1.0e+3
    spec.paramod.z *= u2si
    spec.paramod.h *= u2si
    spec.paramod.rho *= u2si
    spec.paramod.beta *= u2si
    spec.paramod.alpha *= u2si
    spec.paramod.c *= u2si

  if args.winargs is None:
    nwin = 1
    fias = [ colon2lint(args.fidxargs) ]
    cias = [ colon2lint(args.cidxargs) ]
    wwin = [ 1.0 ]
  else:
    nwin = len(args.winargs)
    fias = []
    cias = []
    wwin = []
    for iw in range(nwin):
      winargs = args.winargs[iw].split(',')
      fias.append( colon2lint(winargs[0]) )
      cias.append( colon2lint(winargs[1]) )
      if len(winargs) > 2:
        wwin.append( float(winargs[2]) )
      else:
        wwin.append( 1.0 )

  wwsum = sum(wwin)
  for iw in range(nwin):
    wwin[iw] /= wwsum
    if fias[iw][1] <= 0:
      fias[iw][1] += spec.paramod.fnum
    if cias[iw][1] <= 0:
      cias[iw][1] += spec.paramod.cnum

  wdata = [ data[ fias[iw][0]:fias[iw][1]:fias[iw][2],
    cias[iw][0]:cias[iw][1]:cias[iw][2] ].copy() for iw in range(nwin) ]

  for iw in range(nwin):
    wdata[iw] = tonorm(wdata[iw])

  spec.grtcmod.grtcinitialize()

  #-----------------------------------------------------------------------------

  def tocalc(fname = None):
    if fname:
      print('>> For the file: ' + fname)
      model = np.loadtxt(fname, skiprows = 1)
      spec.paramod.rho   = model[:, 2]
      spec.paramod.beta  = model[:, 3]
      spec.paramod.alpha = model[:, 4]

    misfit = 0.0
    spec.grtcmod.grtcspectrum()
    for iw in range(nwin):
      ker = spec.grtcmod.kermat[ fias[iw][0]:fias[iw][1]:fias[iw][2],
        cias[iw][0]:cias[iw][1]:cias[iw][2] ].copy()
      ker = tonorm(ker)
      diff = ker - wdata[iw]
      misfit += 0.5 * np.sum(diff * diff) / ker.size * wwin[iw]
    return misfit

  #-----------------------------------------------------------------------------

  invfiles = sorted( glob.glob(args.outdir + 'inv*.dat') )
  nitr = len(invfiles)
  mfv = np.zeros((nitr + 1, 2))
  mfv[:, 0] = np.arange(nitr + 1)

  mfv[0, 1] = tocalc()
  for i in range(nitr):
    mfv[i + 1, 1] = tocalc(invfiles[i])
  if not args.Misfile is None:
    np.savetxt(args.Misfile, mfv, fmt = '%5d    %20.10f')

  spec.grtcmod.grtcfinalize()
  spec.paramod.parafinalize()
else:
  mfv = np.loadtxt(args.misfile)

#===============================================================================

maxv = np.amax(mfv[:, 1])
mfv[:, 1] = mfv[:, 1] / maxv

if args.logy:
  plt.semilogy(mfv[:, 0], mfv[:, 1], '-')
  plt.semilogy(mfv[:, 0], mfv[:, 1], 'o', markersize = 3)
else:
  plt.plot(mfv[:, 0], mfv[:, 1], '-')
  plt.plot(mfv[:, 0], mfv[:, 1], 'o', markersize = 3)
plt.gca().xaxis.set_major_locator(tkr.MaxNLocator(integer = True))
plt.xlabel('Iteration number')
plt.ylabel('Normalized misfit')
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
if args.show:
  plt.show()
else:
  plt.savefig(args.savepath + args.outpref + psfx, dpi = args.dpi,
    bbox_inches = 'tight')

