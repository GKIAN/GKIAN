#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Thu 01 Apr 2021 03:43:33 PM CST
#-------------------------------------------------------------------------------

import dispec as spec
from inverse import *
import h5py
import numpy as np
from mpi4py import MPI
from math import ceil
import scipy.optimize as sop
import time, datetime, os
import argparse

#===============================================================================

iopts = { 'iprint': -1,
          'gtol': 1.0e-18, 'ftol': 1.0e-24,
          'maxiter': 999, 'maxfun': 500 }

#===============================================================================

ap = argparse.ArgumentParser(description = 'to inverse S-wave velocity ' +
  'parameters, with dispec.', prefix_chars = '-+')
# for the input files
ap.add_argument('-0', '--refmodel', metavar = 'referenceModel',
  default = 'initial.dat', help = 'the reference model file. ' +
  '[default: %(default)s]')
ap.add_argument('-d', '--datafile', metavar = 'dataFile', default = 'spec.dat',
  help = 'the data input file. [default: %(default)s]')
ag = ap.add_mutually_exclusive_group()
ag.add_argument('-I', '--inml', metavar = 'inputNMLFile', default = 'input.nml',
  help = 'the input NAMELIST file. [default: %(default)s]')
ag.add_argument('-W', '--ifactor', metavar = 'imgFreqFactor', default = None,
  type = float, help = 'the imaginary frequency factor. If not None, ' +
  'read f/c sampling points from the respective input files.')
ap.add_argument('-f', '--ffile', metavar = 'freqFile', default = 'f.txt',
  help = 'the frequency input file. [default: %(default)s]')
ap.add_argument('-c', '--cfile', metavar = 'cFile', default = 'c.txt',
  help = 'the phase velocity input file. [default: %(default)s]')
# for the input interpreter
ap.add_argument('-i', '--fidxargs', metavar = 'iStart[:iEnd[:iStride]]',
  default = '0:0:1', help = 'the arguments for calculating the index of ' +
  'frequency sampling-points in inversion.')
ap.add_argument('-j', '--cidxargs', metavar = 'jStart[:jEnd[:jStride]]',
  default = '0:0:1', help = 'the arguments for calculating the index of ' +
  'phase velocity sampling-points in inversion.')
ap.add_argument('-w', '--winargs', nargs = '+',
  metavar = 'iStart[:iEnd[:iStride]],jStart[:jEnd[:jStride]][,weight]',
  default = None, help = 'the arguments for windowing spectrogram ' +
  'in inversion.')
ap.add_argument('-t', '--datatp', action = 'store_true',
  help = 'transpose the data matrix.')
ap.add_argument('-n', '--nottn', action = 'store_true', help = 'use not ' +
  'the traditional non-exact normalization, but the strict normalization.')
ap.add_argument('-u', '--notsiu', action = 'store_true', help = 'use NOT ' +
  'SI Units in the initial model, but common units, such as km, km/s, g/cm^3.')
# for the inversion controller
ap.add_argument('+i', '--ninit', metavar = 'numInitModel', default = 1,
  type = int, help = 'the number of initial models. [default: %(default)s]')
ap.add_argument('+r', '--rinit', action = 'store_true',
  help = 'randomly generate the different initial models.')
ap.add_argument('-b', '--hwinit', metavar = 'halfWidthInit', type = float,
  default = None, help = "the half width of initial models' range, use m/s " +
  'as unit; if None, it will take min(vs)/4.')
ap.add_argument('+I', '--impinit', metavar = 'fileImpInit', default = None,
  help = 'the name of file for importing initial models. It will ignore ' +
  "'+i', '+r' and '-b'.")
ap.add_argument('-B', '--maxupd', metavar = 'maxUpdate', type = float,
  default = None, help = 'the maximum permissive value for model updating; ' +
  'if None, vs will result from [vs/2, vp].')
ap.add_argument('-E', '--emptype', metavar = 'typeEmpirical', default = 'B',
  help = "the type of empirical relation. B: Brocher's fit; U: User's local " +
  'fit; N: No empirical relation. [default: %(default)s]')
ap.add_argument('-T', '--mistype', metavar = 'typeMisfit', default = 'GCCC',
  help = 'the type of misfit measurement. [default: %(default)s; ' +
  'acceptable: SD2N, GCCC]')
ap.add_argument('-F', '--fwgtarg', metavar = 'fWeightArg', type = float,
  default = None, help = 'the argument for weighting the misfit for ' +
  'different frequencies, only available for SD2N and GCCC. ' +
  '[default: %(default)s]')
ap.add_argument('-r', '--regargs', nargs = 2, metavar = 'regArg', type = float,
  default = [0.0e-9, 0.0e-6], help = 'two regularization arguments for ' +
  '(d2init, smooth). [default: 0.0 0.0]')
# for the output controller
ap.add_argument('-O', '--outpath', metavar = 'outputPath', default = './',
  help = 'save the inversed result(s) to the path. [default: %(default)s]')
ap.add_argument('-v', '--verbosity', action = 'count', default = 0,
  help = 'increase output verbosity.')
args = ap.parse_args()

#===============================================================================

mycomm = MPI.COMM_WORLD
cwsize = mycomm.Get_size()
whoami = mycomm.Get_rank()
isroot = (whoami == 0)

if isroot:
  t0 = time.time()
  ts = datetime.datetime.now().strftime('%a, %d %b %Y, %H:%M:%S %p')
  os.makedirs(args.outpath, exist_ok = True)

if args.datafile[-2:] == 'h5':
  h5id = h5py.File(args.datafile, 'r')
  data = h5id['amp'][:]
  h5id.close()
else:
  data = np.loadtxt(args.datafile)
if args.datatp:
  data = data.T

spec.math.mathinitialize()
if args.ifactor is None:
  if isroot:
    print('>> Reading f/c sampling arguments from the input NAMELIST file.\n',
      flush = True)
  spec.paramod.paragetarguments(args.inml, args.refmodel, deriv = True)
  spec.paramod.parainitialize()
else:
  if isroot:
    print('>> Reading f/c sampling points from the respective input files.\n',
      flush = True)
  spec.paramod.paragetarguments(model = args.refmodel, deriv = True)
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

def colon2lint(string):
  subs = string.split(':')
  ret = [ 0, 0, 1 ]
  for i in range(len(subs)):
    if len(subs[i]) > 0:
      ret[i] = int(subs[i])
  return ret

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

nlyr = spec.paramod.nlayer

wdata = [ data[ fias[iw][0]:fias[iw][1]:fias[iw][2],
  cias[iw][0]:cias[iw][1]:cias[iw][2] ].copy() for iw in range(nwin) ]

xwlev = np.amax(spec.paramod.beta)
if args.maxupd is None:
  lbnd = spec.paramod.beta * ( 1.0 / 2.0 - 1.0 )
  ubnd = spec.paramod.alpha - spec.paramod.beta
  ibnds = [ ( lbnd[i] / xwlev, ubnd[i] / xwlev ) for i in range(nlyr) ]
else:
  ibnds = [ ( - args.maxupd / xwlev, args.maxupd / xwlev) for i in range(nlyr) ]

if args.hwinit is None:
  args.hwinit = np.amin(spec.paramod.beta) / 4.0
args.hwinit /= xwlev

if (data.shape[0] != spec.paramod.fnum) \
  | (data.shape[1] != spec.paramod.cnum):
  raise ValueError('Data shape conflicts with fNum & cNum.')

#===============================================================================

spec.grtcmod.grtcinitialize()

invargs = (spec, args, nwin, fias, cias, wwin, wdata, xwlev, whoami)
if args.mistype == 'SD2N':
  inv = InversionSD2N(invargs)
elif args.mistype == 'GCCC':
  inv = InversionGCCC(invargs)
else:
  inv = Inversion(invargs)
  print('No such type for misfit measurement. [%s]\n' % (args.mistype))

strprint = '>> index: %03d | success: %5r | fun: %10.3e -> %10.3e | ' + \
  'nfev: %3d | njev: %3d | nit: %3d'
rec = { 'i': [], 'f': [], 'x': [], 'p': [] }

def doInv(im, x0):
  rec['i'].append(im)
  # do inversion
  inv.prepare(im, x0)
  mf0 = inv.calc_misfit(x0)
  res = sop.minimize(inv.calc_misfit, x0, method = 'L-BFGS-B',
    jac = inv.calc_gradient, callback = inv.callback,
    bounds = ibnds, options = iopts)
  rec['f'].append( inv.calc_misfit(res.x) )
  rec['x'].append( inv.update_model(res.x, True) )
  # print and output
  inv.save_model(args.outpath + 'invmodel-i%03d.dat' % (im))
  istr = strprint % (im, res.success, mf0, res.fun, res.nfev, res.njev, res.nit)
  if args.verbosity > 0:
    if isroot:
      print('-' * 100 + '\n', flush = True)
    rec['p'].append(istr)
  else:
    print(istr, flush = True)

#===============================================================================

if args.impinit is None:
  ninit = args.ninit
  for im in range(whoami, ninit, cwsize):
    if ninit == 1:
      x0 = np.zeros(nlyr)
    else:
      if args.rinit:
        x0 = np.random.rand(nlyr) * 2.0 - 1.0
      else:
        x0 = ( im / (ninit - 1) * 2.0 - 1.0 ) * np.ones(nlyr)
      x0 *= args.hwinit
    doInv(im, x0)
else:
  minits = np.load(args.impinit)['init']
  ninit = minits.shape[0]
  if minits.shape[1] != nlyr:
    raise IOError("Imcompatible layer number for '+I' and '-0'.")
  for im in range(whoami, ninit, cwsize):
    # inv.spec.paramod.z[:-1] = minits[im, :, 0]
    x0 = (minits[im, :, 1] - inv.xbase) / xwlev
    doInv(im, x0)

if args.verbosity > 0:
  mycomm.Barrier()
  for istr in rec['p']:
    print(istr, flush = True)

#-------------------------------------------------------------------------------

# determine the number of models
nim = len(rec['i'])
gnim = np.array(mycomm.allgather(nim))
# gather all results
sendc = { 'i': gnim, 'f': gnim, 'x': gnim * nlyr }
grec = { 'i': np.empty(ninit, dtype = int), 'f': np.zeros((ninit, 1)),
  'x': np.zeros((ninit, nlyr)) }
for k in sendc.keys():
  mycomm.Gatherv(np.array(rec[k]), recvbuf = (grec[k], sendc[k]), root = 0)
# rearrange the resuls
if isroot:
  sidx = np.argsort(grec['i'])
  grec['f'] = grec['f'][ sidx ]
  grec['x'] = grec['x'][ sidx, : ]

if isroot:
  fout = open(args.outpath + 'allinone.txt', 'wt')
  np.savetxt(fout, np.arange(ninit).reshape((1, ninit)),
    fmt = '    Init-%03d' * ninit)
  np.savetxt(fout, grec['f'].T, fmt = '%12.3e' * ninit)
  np.savetxt(fout, grec['x'].T, fmt = '%12.3f' * ninit)
  fout.close()

if isroot:
  t1 = time.time()
  te = str(datetime.timedelta(seconds = t1 - t0))
  print('\n>> Started at ' + ts + '.')
  print('>> Ended with the elasped time <' + te + '>.\n')

#===============================================================================

spec.grtcmod.grtcfinalize()
spec.paramod.parafinalize()

MPI.Finalize()

# vim:ft=python tw=80 ts=4 sw=2 et ai
