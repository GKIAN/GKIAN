#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Thu 15 Apr 2021 21:49:14 PM CST
#===============================================================================

import numpy as np
import matplotlib.pyplot as plt
import argparse, json
import glob, os, sys

dparfile = 'defpars/model.json'
if not os.path.isfile(dparfile):
  pyhome = os.path.split(sys.argv[0])[0]
  dparfile = pyhome + '/../defpars/model.json'
with open(dparfile) as fin:
  dpars = json.load(fin)

#===============================================================================

ap = argparse.ArgumentParser(description = 'to plot model.')
ap.add_argument('initial', nargs = '?', metavar = 'initModel',
  default = dpars['PA'], help = 'the initial model. [default: %(default)s]')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-1', '--inv', metavar = 'invdModel', default = dpars['-1'],
  help = 'to plot the inversed model. [default: %(default)s]')
ag.add_argument('-R', '--outdir', metavar = 'outputDir', default = dpars['-R'],
  help = 'to plot all models under the inversion-output directory, ' +
  'recursively. [default: %(default)s]')
ap.add_argument('-2', '--real', metavar = 'realModel', default = dpars['-2'],
  help = 'to plot the real model. [default: %(default)s]')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-x', '--idxlist', metavar = 'index', type = int, nargs = '+',
  default = dpars['-x'], help = 'the index list for the inversed models ' +
  'to be plotted. [default: %(default)s]')
ag.add_argument('-X', '--idxarg', metavar = 'iStart:iEnd:iStride',
  default = dpars['-X'], help = 'the arguments for calculating indexes of ' +
  'the inversed models to be plotted. [default: %(default)s]')
ap.add_argument('-b', '--hwinit', metavar = 'halfWidthInit', type = float,
  default = dpars['-b'], help = "the half width of initial models' range, " +
  'use m/s as unit. [default: %(default)s]')
ap.add_argument('-v', '--vlim', metavar = 'vLim', type = float, nargs = 2,
  default = dpars['-v'], help = 'set the velocity limits.')
ap.add_argument('-T', '--ptitle', metavar = 'plotTitle', default = dpars['-T'],
  help = 'set the title for plotting. [default: %(default)s]')
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

if args.idxlist is None:
  if not args.idxarg is None:
    idxarg = [ int(arg) for arg in args.idxarg.split(':') ]
    args.idxlist = list( range(idxarg[0], idxarg[1], idxarg[2]) )

if args.kmunit:
  lu2p = 1.0e-3
else:
  lu2p = 1.0

if args.epsfmt:
  psfx = '.eps'
else:
  psfx = '.png'

#===============================================================================

model = np.loadtxt(args.initial, skiprows = 1)
nlyr = model.shape[0]
z = np.zeros(nlyr + 1)
z[:nlyr] = model[:, 1] * lu2p
z[-1] = 2 * z[-2] - z[-3]
vi = np.zeros(nlyr + 1)
vi[1:] = model[:, 3] * lu2p

d = vi[1:].reshape(nlyr, 1) @ np.ones((1, 3))

if not args.real is None:
  vr = np.zeros(nlyr + 1)
  vr[1:] = np.loadtxt(args.real, skiprows = 1)[:, 3] * lu2p
  d[:, 1] = vr[1:]

#===============================================================================

def toplot(fname = None):
  if not fname is None:
    if fname[-7] == 't':
      ifitr = False
      tid = int(fname[-6:-4])
    else:
      ifitr = True
      tid = int(fname[-7:-4])
    vo = np.zeros(nlyr + 1)
    vo[1:] = np.loadtxt(fname, skiprows = 1)[:, 3] * lu2p
    d[:, 2] = vo[1:]

  plt.figure(figsize = (4.5, 5.5))
  if not args.hwinit is None:
    plt.step(vi - args.hwinit, z, 'y--', where = 'post', label = 'init')
    plt.step(vi + args.hwinit, z, 'y--', where = 'post')
  else:
    plt.step(vi, z, where = 'post', label = 'init')
  if not args.real is None:
    plt.step(vr, z, where = 'post', label = 'real')
  if not fname is None:
    plt.step(vo, z, where = 'post', label = 'invd')
  if args.kmunit:
    plt.xlabel('Velocity (km/s)', fontsize = 15, labelpad = 10.0)
    plt.ylabel('Depth (km)', fontsize = 15)
  else:
    plt.xlabel('Velocity (m/s)', fontsize = 15, labelpad = 10.0)
    plt.ylabel('Depth (m)', fontsize = 15)
  if not args.ptitle is None:
    if len(args.ptitle) > 0:
      plt.title(args.ptitle, fontsize = 18, y = 1.15)
  else:
    if not fname is None:
      if ifitr:
        plt.title('Result for itr = %d' % (tid), fontsize = 18, y = 1.15)
      else:
        plt.title('Result of the thread No.%d' % (tid),
          fontsize = 18, y = 1.15)
  plt.gca().invert_yaxis()
  if not args.vlim is None:
    plt.gca().set_xlim(args.vlim)
  else:
    if not args.hwinit is None:
      vmin = np.amin(d) - args.hwinit
      vmax = np.amax(d) + args.hwinit
    else:
      vmin = np.amin(d)
      vmax = np.amax(d)
    vpad = (vmax - vmin) * 0.1
    plt.gca().set_xlim((vmin - vpad, vmax + vpad))
  plt.gca().set_ylim((max(z), min(z)))
  plt.gca().xaxis.set_ticks_position('top')
  plt.gca().xaxis.set_label_position('top')
  plt.gca().spines['right' ].set_visible(False)
  plt.gca().spines['bottom'].set_visible(False)
  if (not args.real is None) | (not fname is None):
    plt.legend()
  if args.show:
    plt.show()
  else:
    if fname is None:
      plt.savefig(args.savepath + args.outpref + psfx, dpi = args.dpi,
        bbox_inches = 'tight')
    else:
      plt.savefig(args.savepath + args.outpref + '.m%02d%s' % (tid, psfx),
        dpi = args.dpi, bbox_inches = 'tight')
    plt.clf()
    plt.close()

#===============================================================================

if not args.inv is None:
  toplot(args.inv)
elif not args.outdir is None:
  invfiles = sorted( glob.glob(args.outdir + 'inv*.dat') )
  if args.idxlist is None:
    for invfile in invfiles:
      toplot(invfile)
  else:
    for idx in args.idxlist:
      if idx > len(invfiles):
        print('>> No inversed result for the order number <%d>' % (idx))
        continue
      invfile = invfiles[idx - 1]
      toplot(invfile)
else:
  toplot()

