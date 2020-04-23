"""
Hankel transformation using Ogata 2005 method:
H. Ogata, A Numerical Integration Formula Based on the Bessel Functions, Publications of the Research Institute for Mathematical Sciences, vol. 41, no. 4, pp. 949-970, 2005.

Compute Hankel(x, f) = int_0^inf  f(l) J_nv(l*x) ldl

# Usage: 
from hankel_transform import HankelT
HankelT(x, f, phase=0)  # note: nv = abs(phase)

# Note:
the required Nroots and h values depend on the functional form of the integrand. should check convergence by increasing Nroots and decreasing h

Xun Shi 2013
"""

import numpy as np
from scipy import vectorize
from scipy.special import jn_zeros, jn, yn
from numpy import pi,sinh,cosh,tanh


def DEtransform(t):
    return t * tanh(pi / 2. * sinh(t))

def deriv_DEtransform(t):
    res = 1. / cosh(pi / 2. * sinh(t))
    res = pi / 2. * t * cosh(t) * res * res + tanh(pi / 2. * sinh(t))
    return res

def Descrete_Hankel_Transfrom(x, ff, phase=0, Nroots=50, h=0.01):  
# result(x) = int_0^inf  ff(l) J_nv(l*x) ldl
# 1/h: node density; Nroots: # of points to use to approximate the integral
    nv = abs(phase)
    zeros_PI = jn_zeros(nv,Nroots) / pi
    phi_dot = deriv_DEtransform(h * zeros_PI)
    y_k = DEtransform(h*zeros_PI) * pi / h
    w_nv_k = yn(nv, zeros_PI * pi) / jn(nv + 1, zeros_PI * pi)
    J_nv = jn(nv, y_k)
    res = pi / x**2 * w_nv_k * y_k * ff(y_k / x) * J_nv * phi_dot
    return res.sum() * np.sign(phase + 0.1)**phase  # +0.1 to give 1 for phase=0
HankelT = vectorize(Descrete_Hankel_Transfrom)



