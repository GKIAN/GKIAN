
from numpy.distutils.core import Extension, setup

#===============================================================================

CFS = ['-fopenmp']

DFS = []
# DFS.append('DEBUG')
# DFS.append('WITHQ')
DFS.append('IEEE')
DFS.append('COLORPRINT')
DFS.append('DERIV=_BETA') # acceptable value: _BETA, _ALPHA, _RHO

'''
To generate the signature file <test.pyf>, run it on the command line:
  f2py -m dispec -h test.pyf src/constant.F90 src/math.F90 \
    src/parameter.F90 src/grtc.F90
'''

SRCS = ['src/constant.F90', 'src/math.F90', 'src/parameter.F90', 'src/grtc.F90']
#> if need to edit the signature file, after editing, uncomment the follow one,
#> and `make dispec`:
# SRCS.append('test.pyf')

#===============================================================================

for iDF in DFS:
  CFS.append('-D' + iDF)

ext = Extension(
  name = 'dispec',
  sources = SRCS,
  extra_f90_compile_args = CFS,
  extra_link_args = ['-lgomp']
  )

setup(
  name = 'dispec',
  description = 'calculate the multimodal DIspersion SPECtrum for a layered model',
  author = 'Tche LIU',
  author_email = 'tcheliu@mail.ustc.edu.cn',
  version = '0.1.2',
  ext_modules = [ext]
  )
