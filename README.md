# hankel
Python code for fast computation of Hankel transforms using Ogata 2005 method

Ref:
H. Ogata, A Numerical Integration Formula Based on the Bessel Functions, Publications of the Research Institute for Mathematical Sciences, vol. 41, no. 4, pp. 949-970, 2005.

Compute Hankel(x, f(x)) = int_0^inf  f(l) J_nv(l*x) ldl

Usage: 
from hankel_transform import HankelT
HankelT(x,f)

Note:
the required Nroots and h values depend on the functional form of the integrand. should check convergence by increasing Nroots and decreasing h
