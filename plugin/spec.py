#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Sat 17 Apr 2021 20:21:31 PM CST
#===============================================================================

import dispec as spec
from disba import PhaseDispersion
import numpy as np
import h5py
import matplotlib.pyplot as plt
import argparse, glob, os

ap = argparse.ArgumentParser(description = 'to plot data or spectrum.',
  prefix_chars = '-+')
# for basic files
ap.add_argument('-d', '--fdata', metavar = 'dataFile', default = 'spec.dat',
  help = 'the data file. [default: %(default)s]')
ap.add_argument('-0', '--initial', metavar = 'initModel',
  default = 'initial.dat', help = 'the initial model. ' +
  '[default: %(default)s]')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-I', '--inml', metavar = 'inputNMLFile', default = 'input.nml',
  help = 'the input NAMELIST file. [default: %(default)s]')
ag.add_argument('-W', '--ifactor', metavar = 'imgFreqFactor', default = None,
  type = float, help = 'the imaginary frequency factor. If not None, ' +
  'read f/c sampling points from the respective input files.')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-1', '--inv', metavar = 'invdModel', default = None,
  help = 'to plot the inversed model.')
ag.add_argument('-R', '--outdir', metavar = 'outputDir', default = 'output/',
  help = 'to plot all models under the inversion-output directory, ' +
  'recursively. [default: %(default)s]')
ap.add_argument('-2', '--real', metavar = 'realModel', default = 'real.dat',
  help = 'the real model. [default: %(default)s]')
# for frequency and phase velocity
ap.add_argument('-f', '--ffile', metavar = 'freqFile', default = 'f.txt',
  help = 'the frequency input file. [default: %(default)s]')
ap.add_argument('-c', '--cfile', metavar = 'cFile', default = 'c.txt',
  help = 'the phase velocity input file. [default: %(default)s]')
ap.add_argument('-i', '--fidxargs', metavar = 'iStart[:iEnd[:iStride]]',
  default = '0:0:1', help = 'the arguments for calculating the index of ' +
  'frequency sampling-points in inversion.')
ap.add_argument('-j', '--cidxargs', metavar = 'jStart[:jEnd[:jStride]]',
  default = '0:0:1', help = 'the arguments for calculating the index of ' +
  'phase velocity sampling-points in inversion.')
ap.add_argument('-w', '--winargs', nargs = '+',
  metavar = 'iStart[:iEnd[:iStride]],jStart[:jEnd[:jStride]]', default = None,
  help = 'the arguments for windowing spectrogram in inversion.')
ap.add_argument('-r', '--rfistep', metavar = 'rootFreqStep', type = int,
  default = 3, help = 'index step length of frequency sampling ' +
  'for root searching. [default: %(default)s]')
# for data processing
ap.add_argument('-t', '--datatp', action = 'store_true',
  help = 'transpose the data matrix.')
ap.add_argument('-n', '--nottn', action = 'store_true', help = 'use not ' +
  'the traditional non-exact normalization, but the strict normalization.')
ap.add_argument('-u', '--notsiu', action = 'store_true',
  help = 'NOT SI Units used in the initial model, but common units, such as ' +
  'km, km/s, g/cm^3.')
# for supplementary options
ap.add_argument('+f', '--fplot', metavar = 'filePlot', default = 'fplot.dat',
  help = 'the input file for file-item plotting. [default: %(default)s]')
ap.add_argument('+F', '--Fplot', metavar = 'modelFplot', default = None,
  help = 'the model file for file-item plotting; if not None, ' +
  "will re-generate the input file of '+f'.")
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-x', '--idxlist', metavar = 'index', type = int, nargs = '+',
  default = None, help = 'the index list for the inversed models ' +
  'to be plotted.')
ag.add_argument('-X', '--idxargs', metavar = 'iStart:iEnd:iStride',
  default = None, help = 'the arguments for calculating indexes of ' +
  'the inversed models to be plotted.')
# for data selecting to plot
ap.add_argument('+b', '--bgitem', metavar = 'itemName', default = None,
  help = 'the item for window background pcolor. [acceptable: data, init, ' +
  'real, inv, file]')
ap.add_argument('+c', '--fgcont', metavar = 'itemName', default = None,
  help = 'the item for window foreground contour. [acceptable: data, init, ' +
  'real, inv, file]')
ap.add_argument('+p', '--fgplot', metavar = 'itemName', default = None,
  help = 'the item for window foreground plot. [acceptable: data, init, ' +
  'real, inv, file]')
# for plotting options
ap.add_argument('+B', '--bgmain', action = 'store_true',
  help = 'to pcolor the main background, use itemName same as bgitem.')
ap.add_argument('-L', '--clev', metavar = 'level', type = float, nargs = '+',
  default = [0.25, 0.75], help = 'the level for contour plotting. ' +
  '[default: 0.25 0.75]')
ap.add_argument('-T', '--ptitle', metavar = 'plotTitle', default = None,
  help = 'set the title for plotting.')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-s', '--savepath', metavar = 'saveFigPath', default = './',
  help = 'save figure(s) to the path. [default: %(default)s]')
ag.add_argument('-S', '--show', action = 'store_true',
  help = 'show figure(s) immediately.')
ap.add_argument('-O', '--outpref', metavar = 'outFigPrefix',
  default = 'Ex4I_DataSpec', help = 'save figure as a file with ' +
  'the name prefix. [default: %(default)s]')
ap.add_argument('-p', '--dpi', metavar = 'dpiValue', type = int, default = 150,
  help = 'specify the value of dpi for saving figure. [default: %(default)s]')
ap.add_argument('-k', '--kmunit', action = 'store_true',
  help = 'use km as the length unit.')
ap.add_argument('-e', '--epsfmt', action = 'store_true',
  help = 'save figure(s) with EPS format; otherwise, with PNG format.')
args = ap.parse_args()

if args.idxlist is None:
  if not args.idxargs is None:
    idxargs = [ int(ia) for ia in args.idxargs.split(':') ]
    args.idxlist = list( range(idxargs[0], idxargs[1], idxargs[2]) )

if args.kmunit:
  lu2p = 1.0e-3
else:
  lu2p = 1.0

if args.epsfmt:
  psfx = '.eps'
else:
  psfx = '.png'

#===============================================================================

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
else:
  u2si = 1.0
spec.paramod.z *= u2si
spec.paramod.h *= u2si
spec.paramod.c *= u2si

def colon2lint(string):
  subs = string.split(':')
  ret = [ 0, 0, 1 ]
  for i in range(len(subs)):
    ret[i] = int(subs[i])
  return ret

if args.winargs is None:
  nwin = 1
  fias = [ colon2lint(args.fidxargs) ]
  cias = [ colon2lint(args.cidxargs) ]
else:
  nwin = len(args.winargs)
  fias = []
  cias = []
  for iw in range(nwin):
    winargs = args.winargs[iw].split(',')
    fias.append( colon2lint(winargs[0]) )
    cias.append( colon2lint(winargs[1]) )

for iw in range(nwin):
  if fias[iw][1] <= 0:
    fias[iw][1] += spec.paramod.fnum
  if cias[iw][1] <= 0:
    cias[iw][1] += spec.paramod.cnum

fs = spec.paramod.f.copy()
cs = spec.paramod.c.copy()

wfs = [ fs[ fias[iw][0]:fias[iw][1]:fias[iw][2] ] for iw in range(nwin) ]
wcs = [ cs[ cias[iw][0]:cias[iw][1]:cias[iw][2] ] for iw in range(nwin) ]

model = np.zeros((spec.paramod.nlayer, 5))

itfile = dict({ 'data': args.fdata, 'init': args.initial, 'real': args.real,
  'file': args.fplot, 'inv': args.inv })
itvwin = dict({ 'data': None, 'init': None, 'real': None, 'file': None,
  'inv': None })
if args.bgmain:
  Bmat = np.zeros((spec.paramod.fnum, spec.paramod.cnum))

for item in [args.bgitem, args.fgcont, args.fgplot]:
  if (not item is None) & (not item in itfile.keys()):
    raise ValueError("Not acceptable value <" + item + "> for 'itemName'.")
if not args.fgplot is None:
  if (args.fgplot == args.bgitem) | (args.fgplot == args.fgcont):
    raise ValueError("The option '+p' conflict with '+b' or '+c'.")

spec.grtcmod.grtcinitialize()

#===============================================================================

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

def model2data(fname, srs = 0, isspec = True, ufctr = 1.0):
  model[:, :] = np.loadtxt(fname, skiprows = srs)[:, :5]
  if isspec:
    spec.paramod.rho   = model[:, 2] * ufctr
    spec.paramod.beta  = model[:, 3] * ufctr
    spec.paramod.alpha = model[:, 4] * ufctr
    spec.grtcmod.grtcspectrum()
    var = spec.grtcmod.kermat.copy()
    var = tonorm(var, args.nottn)
  else:
    model[:, 1:] *= ufctr * 1.0e-3
    model[:, 1] = np.append( np.diff(model[:, 1]), [0,] )
    pD = PhaseDispersion(*model[:, [1, 4, 3, 2]].T, algorithm = 'dunkin')
    nmode = 10
    rfcs = [ [], [] ]
    for iw in range(nwin):
      for m in range(nmode):
          iwfs = wfs[iw][::args.rfistep]
          cpR = pD(1.0 / iwfs[::-1], m, wave = 'rayleigh')
          fR = 1.0 / cpR.period[::-1]
          cR = cpR.velocity[::-1] * 1.0e+3
          inWin = (cR >= wcs[iw][0]) & (cR <= wcs[iw][-1])
          rfcs[0] = np.append(rfcs[0], fR[inWin])
          rfcs[1] = np.append(rfcs[1], cR[inWin])
    var = np.array(rfcs).T
  return var

def tocalc(item, fname = None):
  srs = 1
  if not fname is None:
    itfile[item] = fname
  if (args.bgitem == item) | (args.fgcont == item):
    isspec = True
  elif args.fgplot == item:
    isspec = False
  else:
    return
  # to calculate
  if item in ['init', 'real']:
    res = model2data(itfile[item], srs, isspec, u2si)
  elif item == 'inv':
    res = model2data(itfile[item], srs, isspec)
  elif item == 'data':
    if args.fdata[-2:] == 'h5':
      h5id = h5py.File(itfile[item], 'r')
      res = h5id['amp'][:]
      h5id.close()
    else:
      res = np.loadtxt(itfile[item])
    if args.datatp:
      res = res.T
  elif item == 'file':
    if args.Fplot is None:
      res = np.loadtxt(itfile[item])
    else:
      res = model2data(args.Fplot, 1, isspec, 1.0)
      np.savetxt(itfile[item], res)
  if args.bgmain & (args.bgitem == item):
    Bmat[:, :] = tonorm(res, args.nottn)
  # to window
  if isspec:
    itvwin[item] = [ res[ fias[iw][0]:fias[iw][1]:fias[iw][2],
      cias[iw][0]:cias[iw][1]:cias[iw][2] ].copy() for iw in range(nwin) ]
    for iw in range(nwin):
      itvwin[item][iw] = tonorm(itvwin[item][iw], args.nottn)
  else:
    itvwin[item] = res

def toplot(finv = None):
  if not finv is None:
    if finv[-8] == 'i':
      ifitr = False
    else:
      ifitr = True
    tid = int(finv[-7:-4])
    tocalc('inv', finv)
  # main plot
  if args.bgmain:
    plt.pcolor(fs, cs * lu2p, Bmat.T, cmap = 'jet', shading = 'auto',
      rasterized = True, alpha = 0.5, snap = True)
  if not args.bgitem is None:
    for iw in range(nwin):
      plt.pcolormesh(wfs[iw], wcs[iw] * lu2p, itvwin[args.bgitem][iw].T,
        cmap = 'jet', shading = 'auto', rasterized = True)
        # cmap = 'jet', shading = 'gouraud', rasterized = True)
    plt.colorbar()
  if not args.fgcont is None:
    for iw in range(nwin):
      plt.contour(wfs[iw], wcs[iw] * lu2p, itvwin[args.fgcont][iw].T,
        levels = args.clev)
    if args.bgitem is None:
      plt.colorbar()
  if not args.fgplot is None:
    if args.bgitem is None:
      lstyle = 'k.'
    else:
      lstyle = 'w.'
    plt.plot(itvwin[args.fgplot][:, 0], itvwin[args.fgplot][:, 1] * lu2p,
      lstyle, markersize = 2)
  # set axies
  if not args.bgmain:
    plt.gca().set_xlim((fs[0], fs[-1]))
    plt.gca().set_ylim((cs[0] * lu2p, cs[-1] * lu2p))
  plt.xlabel('Frequency (Hz)')
  if args.kmunit:
    plt.ylabel('Phase velocity (km/s)')
  else:
    plt.ylabel('Phase velocity (m/s)')
  if not args.ptitle is None:
    plt.title(args.ptitle)
  elif not finv is None:
    if ifitr:
      plt.title('Result for itr = %d' % (tid))
    else:
      plt.title('Result of the initial-model No.%d' % (tid))
  # display
  if args.show:
    plt.show()
  else:
    if finv is None:
      plt.savefig(args.savepath + args.outpref + psfx, dpi = args.dpi,
        bbox_inches = 'tight')
    else:
      plt.savefig(args.savepath + args.outpref + '.m%02d%s' % (tid, psfx),
        dpi = args.dpi, bbox_inches = 'tight')
    plt.clf()

#===============================================================================

for item in list(itvwin.keys())[:-1]:
  tocalc(item)

if (args.bgitem == 'inv') | (args.fgcont == 'inv') | (args.fgplot == 'inv'):
  if not args.inv is None:
    toplot(args.inv)
  else:
    invfiles = sorted( glob.glob(args.outdir + 'inv*.dat') )
    if args.idxlist is None:
      for invfile in invfiles:
        print('>> For the file: ' + invfile)
        toplot(invfile)
    else:
      for idx in args.idxlist:
        if idx > len(invfiles):
          print('>> No inversed result for the order number <%d>' % (idx))
          continue
        invfile = invfiles[idx - 1]
        print('>> For the file: ' + invfile)
        toplot(invfile)
else:
  toplot()

#===============================================================================

spec.grtcmod.grtcfinalize()
spec.paramod.parafinalize()
