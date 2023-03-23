#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Sun 26 Jun 2022 09:46:18 AM CST
#-------------------------------------------------------------------------------

from disba import PhaseDispersion
import numpy as np
import matplotlib.pyplot as plt
import argparse, glob, os

ap = argparse.ArgumentParser(description = 'to plot data or curve.',
  prefix_chars ='-+')
ap.add_argument('-d', '--fdata', metavar = 'dataFile', default = None,
  help = 'the data file.')
ap.add_argument('-0', '--initial', metavar = 'initModel', default = None,
  help = 'the initial model.')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-1', '--inv', metavar = 'invdModel', default = None,
  help = 'to plot the inversed model.')
ag.add_argument('-A', '--aiofile', metavar = 'allInOneFile', default = None,
  help = 'to plot some inversed models for different threads.')
ap.add_argument('-2', '--real', metavar = 'realModel', default = None,
  help = 'to plot the real model.')
ap.add_argument('-a', '--avrg', metavar = 'avrgModel', default = None,
  help = 'to plot the averaged model.')
ap.add_argument('-R', '--outdir', metavar = 'outputDir', default = 'output/',
  help = 'the inversion-output directory. [default: %(default)s]')
ap.add_argument('-f', '--ffile', metavar = 'freqFile', default = None,
  help = 'the frequency input file.')
ap.add_argument('-c', '--clim', metavar = 'cLim', type = float, nargs = 2,
  default = None, help = 'set the phase-velocity limits.')
ap.add_argument('-m', '--nmode', metavar = 'modeNum', type = int,
  default = None, help = 'the number of mode to search. If has fdata, ' +
  'take the maximum value of modes in fdata.')
ap.add_argument('-P', '--percplot', metavar = 'percentPlot', type = float,
  default = 100, help = 'use the best percentage of all models to plot for ' +
  'allinone. [default: %(default)s]')
ap.add_argument('-i', '--fidxargs', metavar = 'iStart[:iEnd[:iStride]]',
  default = '0:0:1', help = 'the arguments for calculating the index of ' +
  'frequency sampling-points in inversion.')
ap.add_argument('-l', '--pline', action = 'store_true',
  help = 'use connected-line, NOT scattered-markers to plot.')
ap.add_argument('-u', '--notsiu', action = 'store_true',
  help = 'NOT SI Units used in forward modelling, but common units, such as ' +
  'km, km/s, g/cm^3.')
ap.add_argument('-T', '--ptitle', metavar = 'plotTitle', default = None,
  help = 'set the title for plotting.')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-s', '--savepath', metavar = 'saveFigPath', default = './',
  help = 'save figure(s) to the path. [default: %(default)s]')
ag.add_argument('-S', '--show', action = 'store_true',
  help = 'show figure(s) immediately.')
ap.add_argument('-O', '--outpref', metavar = 'outFigPrefix',
  default = 'Ex4I_DataCurve', help = 'save figure as a file with ' +
  'the name prefix. [default: %(default)s]')
ap.add_argument('-p', '--dpi', metavar = 'dpiValue', type = int, default = 150,
  help = 'specify the value of dpi for saving figure. [default: %(default)s]')
ap.add_argument('-k', '--kmunit', action = 'store_true',
  help = 'use km as the length unit.')
ap.add_argument('-e', '--epsfmt', action = 'store_true',
  help = 'save figure(s) with EPS format; otherwise, with PNG format.')
args = ap.parse_args()

if args.notsiu:
  u2km = 1.0
else:
  u2km = 1.0e-3

if args.kmunit:
  lu2p = 1.0
else:
  lu2p = 1.0e+3

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

if not args.fdata is None:
  data = np.loadtxt(args.fdata)
  data[:, 1] *= u2km
  fs = data[:, 0]
  if data.shape[1] > 2:
    args.nmode = 1 + np.amax(data[:, 2]).astype(int)
    dfcs = [ [] for m in range(args.nmode) ]
    for f, c, x in data:
      m = x.astype(int)
      dfcs[m].append([ f, c ])
    for m in range(args.nmode):
      dfcs[m] = np.array(dfcs[m])
  else:
    dfcs = [ data ]

fias = colon2lint(args.fidxargs)
if not args.ffile is None:
  fs = np.loadtxt(args.ffile)
  if fias[1] <= 0:
    fias[1] += len(fs)
  fs = fs[ fias[0]:fias[1]:fias[2] ]

if not 'fs' in vars():
  raise SyntaxError("'-d' or/and '-f' MUST be specified.")

def model2data(fname, srs = 0, ufactor = u2km):
  model = np.loadtxt(fname, skiprows = srs)[:, [1, 4, 3, 2]] * ufactor
  model[:, 0] = np.append( np.diff(model[:, 0]), [0,] )
  pD = PhaseDispersion(*model.T, algorithm = 'dunkin')
  if not args.nmode  is None:
    vfcs = [ [] for m in range(args.nmode) ]
    for m in range(args.nmode):
      cpR = pD(1.0 / fs[::-1], m, wave = 'rayleigh')
      fR = 1.0 / cpR.period[::-1]
      cR = cpR.velocity[::-1]
      vfcs[m] = np.vstack((fR, cR)).T
    return vfcs
  else:
    nmode = 10
    vfcs = [ [], [] ]
    for m in range(nmode):
      cpR = pD(1.0 / fs[::-1], m, wave = 'rayleigh')
      fR = 1.0 / cpR.period[::-1]
      cR = cpR.velocity[::-1]
      vfcs[0] = np.append(vfcs[0], fR)
      vfcs[1] = np.append(vfcs[1], cR)
    var = [ np.array(vfcs).T ]
    return var

if args.pline:
  cresp = 'color'
  comm = { 'linestyle': '-' }
else:
  cresp = 'markerfacecolor'
  comm = { 'linestyle': 'None', 'marker': '.', 'markersize': 5,
    'markeredgecolor': 'None' }

resp = {
  'data': { cresp: 'turquoise' },
  'init': { cresp: '#BFBF00' },
  'invd': { cresp: 'black' },
  'real': { cresp: 'teal' },
  'avrg': { cresp: 'red' }
  }

def toplot(fcs, who, label = True, weight = 1.0):
  nmode = len(fcs)
  for m in range(nmode):
    if fcs[m].shape[0] == 0:
      continue
    if not args.fdata is None:
      fi0 = np.argmin(np.abs(fcs[m][:, 0] - dfcs[m][ 0, 0]))
      fi1 = np.argmin(np.abs(fcs[m][:, 0] - dfcs[m][-1, 0]))
    else:
      fi0 =  0
      fi1 = -1
    if label:
      plt.plot(fcs[m][fi0:fi1, 0], fcs[m][fi0:fi1, 1] * lu2p, alpha = weight,
        **comm, **resp[who], label = who)
      label = False
    else:
      plt.plot(fcs[m][fi0:fi1, 0], fcs[m][fi0:fi1, 1] * lu2p, alpha = weight,
        **comm, **resp[who])

#===============================================================================

if not args.fdata is None:
  toplot(dfcs, 'data')

if not args.initial is None:
  vfcs = model2data(args.initial, 1)
  toplot(vfcs, 'init')

if not args.inv is None:
  vfcs = model2data(args.inv, 1)
  toplot(vfcs, 'invd')
elif not args.aiofile is None:
  mfvs = np.loadtxt(args.aiofile, ndmin = 2, skiprows = 1)[0, :]
  npick = int( len(mfvs) * args.percplot / 100 )
  if npick < 1:
    raise ValueError('Too small best-percentage <%g> to plot.' % (args.percplot))
  ipick = np.argsort(mfvs)[:npick]
  mpick = mfvs[ipick]
  mpick -= np.amin(mpick)
  if np.all(mpick == 0.0):
    wpick = np.ones_like(mpick)
  else:
    mpick /= np.amax(mpick)
    wpick = np.exp( - mpick)
  wpick /= np.sum(wpick)
  label = True
  for i in range(npick):
    ifile = args.outdir + 'invmodel-i%03d.dat' % (ipick[i])
    vfcs = model2data(ifile, 1)
    toplot(vfcs, 'invd', label, wpick[i])
    if label:
      label = False

if not args.real is None:
  vfcs = model2data(args.real, 1)
  toplot(vfcs, 'real')

if not args.avrg is None:
  vfcs = model2data(args.avrg, 1)
  toplot(vfcs, 'avrg')

plt.legend()
if not args.clim is None:
  plt.gca().set_ylim((args.clim))

plt.xlabel('Frequency (Hz)')
if args.kmunit:
  plt.ylabel('Phase velocity (km/s)')
else:
  plt.ylabel('Phase velocity (m/s)')
if not args.ptitle is None:
  plt.title(args.ptitle)

if args.show:
  plt.show()
else:
  plt.savefig(args.savepath + args.outpref + psfx, dpi = args.dpi,
    bbox_inches = 'tight')

# vim:ft=python tw=80 ts=4 sw=2 et ai
