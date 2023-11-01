#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Tue 31 Oct 2023 08:39:49 PM CST
#-------------------------------------------------------------------------------

import json, sys

#===============================================================================

if 1:
  # for topy/pyinv:
  outfile = '../defpars/pyinv.json'
  pars = {
      '-0': 'initial.dat',
      '-d': 'spec.dat',
      '-I': 'input.nml',
      '-W': None,
      '-f': 'f.txt',
      '-c': 'c.txt',
      '-i': '0:0:1',
      '-j': '0:0:1',
      '-w': None,
      '+i': 1,
      '-b': None,
      '+I': None,
      '-B': None,
      '-E': 'B',
      '-T': 'GCCC',
      '-F': None,
      '-r': [0.0, 0.0],
      '-O': './',
      '-v': 0,
    }

if 0:
  # for topy/callc:
  outfile = '../defpars/callc.json'
  pars = {
      '-0': 'initial.dat',
      '-d': 'spec.dat',
      '-I': 'input.nml',
      '-W': None,
      '-f': 'f.txt',
      '-c': 'c.txt',
      '-i': '0:0:1',
      '-j': '0:0:1',
      '-w': None,
      '-B': None,
      '-E': 'B',
      '-T': 'GCCC',
      '-F': None,
      '-r': [0.0, 0.0],
      '-l': '@-7:0:0',
      '-L': 'S',
      '-O': './',
      '-v': 0,
    }

#-------------------------------------------------------------------------------

if 0:
  # for plugin/alltone:
  outfile = '../defpars/alltone.json'
  pars = {
      'PA': 'output/allione.txt',
      '-0': 'initial.dat',
      '-2': None,
      '-E': 'B',
      '-v': None,
      '-P': 50,
      '-a': None,
      '-b': None,
      '-V': None,
      '-s': './',
      '-O': 'Ex4I_Result',
      '-p': 150,
    }

if 0:
  # for plugin/curve:
  outfile = '../defpars/curve.json'
  pars = {
      '-d': None,
      '-0': None,
      '-1': None,
      '-A': None,
      '-2': None,
      '-a': None,
      '-R': 'output/',
      '-f': None,
      '-c': None,
      '-m': None,
      '-P': 100,
      '-i': '0:0:1',
      '-T': None,
      '-s': './',
      '-O': 'Ex4I_DataCurve',
      '-p': 150,
    }

if 0:
  # for plugin/misfit:
  outfile = '../defpars/misfit.json'
  pars = {
      '-m': None,
      'EP': 'spec.dat',
      '-0': 'initial.dat',
      '-I': 'input.nml',
      '-W': None,
      '-f': 'f.txt',
      '-c': 'c.txt',
      '-R': 'output/',
      '-i': '0:0:1',
      '-j': '0:0:1',
      '-w': None,
      '-M': None,
      '-s': './',
      '-O': 'Ex4I_Misfit',
      '-p': 150,
    }

if 0:
  # for plugin/model:
  outfile = '../defpars/model.json'
  pars = {
      'PA': 'initial.dat',
      '-1': None,
      '-R': None,
      '-2': None,
      '-x': None,
      '-X': None,
      '-b': None,
      '-v': None,
      '-T': None,
      '-s': './',
      '-O': 'Ex4I_Model',
      '-p': 150,
    }

if 0:
  # for plugin/senker:
  outfile = '../defpars/senker.json'
  pars = {
      '+k': None,
      'EP': 'model.dat',
      '-I': 'input.nml',
      '-W': None,
      '-f': 'f.txt',
      '-c': 'c.txt',
      '-i': '0:0:1',
      '-j': '0:0:1',
      '-w': None,
      '+K': None,
      '-E': None,
      '-s': './',
      '-O': 'Ex4I_Kernel',
      '-p': 150,
    }

if 0:
  # for plugin/spec:
  outfile = '../defpars/spec.json'
  pars = {
      '-d': 'spec.dat',
      '-0': 'initial.dat',
      '-I': 'input.nml',
      '-W': None,
      '-1': None,
      '-R': 'output/',
      '-2': 'real.dat',
      '-f': 'f.txt',
      '-c': 'c.txt',
      '-i': '0:0:1',
      '-j': '0:0:1',
      '-w': None,
      '-r': 3,
      '+f': 'fplot.dat',
      '+F': None,
      '-x': None,
      '-X': None,
      '+b': None,
      '+c': None,
      '+p': None,
      '-L': [0.25, 0.75],
      '-T': None,
      '-s': './',
      '-O': 'Ex4I_DataSpec',
      '-p': 150
    }

#===============================================================================

with open(outfile, 'w') as fout:
  json.dump(pars, fout, indent = 2)
  fout.write('\n')

# vim:ft=python tw=80 ts=4 sw=2 et ai
