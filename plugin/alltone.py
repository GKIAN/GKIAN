#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Thu 29 Apr 2021 09:30:12 AM CST
#-------------------------------------------------------------------------------

import empirical as emp
import numpy as np
import scipy.stats as sst
import matplotlib.pyplot as plt
import argparse

ap = argparse.ArgumentParser(description = 'to calculate and plot ' +
  'from all to one.')
ap.add_argument('aiofile', metavar = 'allInOneFile', nargs = '?',
  default = 'output/allinone.txt', help = 'the AllInOne file. ' +
  '[default: %(default)s]')
ap.add_argument('-0', '--initial', metavar = 'initModel',
  default = 'initial.dat', help = 'the initial model. [default: %(default)s]')
ap.add_argument('-2', '--real', metavar = 'realModel', default = None,
  help = 'to plot the real model.')
ap.add_argument('-E', '--emptype', metavar = 'typeEmpirical', default = 'B',
  help = "the type of empirical relation. B: Brocher's fit, U: User's local " +
  'fit, N: No empirical relation. [default: %(default)s]')
ap.add_argument('-v', '--vlim', metavar = 'vLim', nargs = 2, type = float,
  default = None, help = 'the velocity limits (vmin, vmax).')
ap.add_argument('-P', '--percplot', metavar = 'percentPlot', type = float,
  default = 50, help = 'use the best percentage of all models to plot. ' +
  '[default: %(default)s]')
ap.add_argument('-a', '--avrgfile', metavar = 'averageFile', default = None,
  help = 'the output filename of average model; if None, not output.')
ap.add_argument('-b', '--hwinit', metavar = 'halfWidthInit', type = float,
  default = None, help = "the half width of initial models' range, use m/s " +
  'as unit.')
ap.add_argument('-K', '--usekde', action = 'store_true',
  help = 'use KDE statistical analysis.')
ap.add_argument('-V', '--vmaxkde', metavar = 'vMaxKDE', type = float,
  default = None, help = 'the vmax value of pcolor for KDE.')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-s', '--savepath', metavar = 'saveFigPath', default = './',
  help = 'save figure(s) to the path. [default: %(default)s]')
ag.add_argument('-S', '--show', action = 'store_true',
  help = 'show figure(s) immediately.')
ap.add_argument('-O', '--outpref', metavar = 'outFigPrefix',
  default = 'Ex4I_Result', help = 'save figure as a file with ' +
  'the name prefix. [default: %(default)s]')
ap.add_argument('-p', '--dpi', metavar = 'dpiValue', type = int, default = 150,
  help = 'specify the value of dpi for saving figure. [default: %(default)s]')
ap.add_argument('-k', '--kmunit', action = 'store_true',
  help = 'use km as the length unit.')
ap.add_argument('-e', '--epsfmt', action = 'store_true',
  help = 'save figure(s) with EPS format; otherwise, with PNG format.')
args = ap.parse_args()

ngrid = 500

if args.kmunit:
  lu2p = 1.0e-3
else:
  lu2p = 1.0

if args.epsfmt:
  psfx = '.eps'
else:
  psfx = '.png'

#===============================================================================

def calc_pdf(xorg, xgrd, wt = None):
  if len(xorg) == 1:
    ret = np.zeros_like(xgrd)
    idx = np.argmin( np.abs(xorg[0] - xgrd) )
    ret[idx] = 1.0
  else:
    if np.allclose(xorg, xorg[0]):
      x = xorg * ( 1.0 + (np.random.rand(len(xorg)) - 0.5) * 0.10 )
    else:
      x = xorg
    kde = sst.gaussian_kde(x, weights = wt)
    ret = kde.evaluate(xgrd)
  return ret

def update_model(vs):
  model = np.loadtxt(args.initial, skiprows = 1)[:, :5]
  beta = vs / lu2p
  if args.emptype == 'B':
    model[:, 2:5] = emp.brocherFit(beta)
  elif args.emptype == 'U':
    model[:, 2:5] = emp.userFit(beta)
  else:
    model[:, 3] = beta
  return model

#===============================================================================

res = np.loadtxt(args.aiofile, ndmin = 2, skiprows = 1)
mfvs = res[0, :]
vres = res[1:, :] * lu2p
nlyr, nres = vres.shape
npick = int( nres * args.percplot / 100 )
if npick < 1:
  raise ValueError('Too small best-percentage <%g> to plot.' % (args.percplot))

if not args.hwinit is None:
  vinit = np.zeros(nlyr + 1)
  vinit[1:] = np.loadtxt(args.initial, skiprows = 1)[:, 3] * lu2p
  args.hwinit *= lu2p

if not args.real is None:
  vreal = np.zeros(nlyr + 1)
  vreal[1:] = np.loadtxt(args.real, skiprows = 1)[:, 3] * lu2p

ipick = np.argsort(mfvs)[:npick]
mpick = mfvs[ipick]
vpick = vres[:, ipick]

if args.vlim is None:
  vmin = np.amin(vpick)
  vmax = np.amax(vpick)
  if not args.hwinit is None:
    vmin = min(vmin, np.amin(vinit[1:] - args.hwinit))
    vmax = max(vmax, np.amax(vinit[1:] + args.hwinit))
  if not args.real is None:
    vmin = min(vmin, np.amin(vreal[1:]))
    vmax = max(vmax, np.amax(vreal[1:]))
  vpad = (vmax - vmin) * 0.1
  args.vlim = [ vmin - vpad, vmax + vpad ]

mpick -= np.amin(mpick)
if np.all(mpick == 0.0):
  wpick = np.ones_like(mpick)
else:
  mpick /= np.amax(mpick)
  wpick = np.exp( - mpick)
wpick /= np.sum(wpick)

if args.usekde:
  vgrid = np.linspace(args.vlim[0], args.vlim[1], ngrid)
  vpdf = np.zeros((nlyr, ngrid))
  for i in range(nlyr):
    vpdf[i, :] = calc_pdf(vpick[i, :], vgrid, wpick)
  vgrid = np.append(vgrid, args.vlim[1])

z = np.loadtxt(args.initial, skiprows = 1)[:, 1] * lu2p
z = np.append(z, 2 * z[-1] - z[-2])
vavrg = np.zeros(nlyr + 1)
vavrg[1:] = vpick @ wpick
if not args.avrgfile is None:
  mavrg = update_model(vavrg[1:])
  fout = open(args.avrgfile, 'wt')
  fout.write('Layer-No.   Depth(m)   Densities(kg/m3)  S-Velocity(m/s)  ' +
    'P-Velocity(m/s)\n')
  np.savetxt(fout, mavrg, fmt = '%5d      %10.2f   %14.3f    ' +
    '%13.3f    %13.3f')
  fout.close()

#===============================================================================

plt.figure(figsize = (4.5, 5.5))

if not args.hwinit is None:
  plt.step(vinit - args.hwinit, z, 'y--', where = 'post', label = 'init')
  plt.step(vinit + args.hwinit, z, 'y--', where = 'post')

if args.usekde:
  if args.vmaxkde is None:
    plt.pcolormesh(vgrid, z, vpdf, cmap = 'Greys', shading = 'auto')
    plt.colorbar()
  else:
    plt.pcolormesh(vgrid, z, vpdf, cmap = 'Greys', shading = 'auto',
      vmax = args.vmaxkde)
else:
  wmin = np.amin(wpick)
  wmax = np.amax(wpick)
  if wmin != wmax:
    wpick = (wpick - wmin) / (wmax - wmin)
    wpick = wpick * 0.8 + 0.1
  vplot = np.zeros(nlyr + 1)
  vplot[1:] = vpick[:, 0]
  plt.step(vplot, z, 'k-', where = 'post', alpha = wpick[0], label = 'invd')
  for i in range(1, npick):
    vplot[1:] = vpick[:, i]
    plt.step(vplot, z, 'k-', where = 'post', alpha = wpick[i])

if npick > 1:
  plt.step(vavrg, z, 'red', linestyle = 'dotted', where = 'post',
    label = 'avrg')
if not args.real is None:
  plt.step(vreal, z, 'teal', linestyle = '-.', where = 'post',
    label = 'real')

plt.legend()
if args.kmunit:
  plt.xlabel('Velocity (km/s)', fontsize = 15, labelpad = 10.0)
  plt.ylabel('Depth (km)', fontsize = 15)
else:
  plt.xlabel('Velocity (m/s)', fontsize = 15, labelpad = 10.0)
  plt.ylabel('Depth (m)', fontsize = 15)
plt.gca().invert_yaxis()
plt.gca().set_xlim(args.vlim)
plt.gca().set_ylim((max(z), min(z)))
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.gca().spines['right' ].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
if args.show:
  plt.show()
else:
  plt.savefig(args.savepath + args.outpref + psfx, dpi = args.dpi,
    bbox_inches = 'tight')

# vim:ft=python tw=80 ts=4 sw=2 et ai
