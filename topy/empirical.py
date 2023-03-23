
#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Tue 23 Aug 2022 04:50:57 PM CST
#-------------------------------------------------------------------------------

import numpy as np

def brocherFit(vs):
  nlyr = len(vs)
  rba = np.zeros((nlyr, 3))
  for i in range(nlyr):
    beta = vs[i]
    alpha = 0.9409e+03 + 2.0947e+00 * beta - 0.8206e-03 * beta * beta \
      + 0.2683e-06 * (beta ** 3) - 0.0251e-09 * (beta ** 4)
    rho = 1.6612e+00 * alpha - 0.4721e-03 * alpha * alpha \
      + 0.0671e-06 * (alpha ** 3) - 0.0043e-09 * (alpha ** 4) \
      + 0.000106e-12 * (alpha ** 5)
    rba[i, :] = rho, beta, alpha
  return rba

def userFit(vs):
  # Now with the Gardner's rule
  nlyr = len(vs)
  rba = np.zeros((nlyr, 3))
  for i in range(nlyr):
    beta = vs[i]
    alpha = 1.73 * beta
    rho = 1.74e+3 * (alpha / 1.0e+3) ** 0.25
    rba[i, :] = rho, beta, alpha
  return rba

# vim:ft=python tw=80 ts=4 sw=2 et ai
