# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 10:02:16 2022

@author: Pena
"""

from numpy import radians
from sympy import sin, cos, pi, Abs, Heaviside, Integral
import numpy as np
import sympy as smp
import scipy as sp
import matplotlib.pyplot as plt

#from sympy.physics.vector import *
from sympy.vector import CoordSys3D, dot
T = CoordSys3D('T')
r = 1
alphaDeg = 23.5
w = 2*pi/365
W = 732*pi/365
plot = 0
plotArrows = 0

alpha = radians(alphaDeg)
sinAlpha = sin(alpha)
cosAlpha = cos(alpha)

#theta, phi, t, w, W, alpha = smp.symbols('theta, phi, t, omega, Omega, alpha')
theta, phi, t, alpha = smp.symbols('theta, phi, t, alpha', real=True)

P = sin(theta)*cos(phi)*T.i + sin(theta)*sin(phi)*T.j + cos(theta)*T.k
S = (cos(alpha)*sin(W*t)*sin(w*t) + cos(W*t)*cos(w*t))*T.i + \
    (cos(alpha)*cos(W*t)*sin(w*t) - sin(W*t)*cos(w*t))*T.j + \
    (sin(alpha)*sin(w*t))*T.k

f = abs(dot(S, P))*sin(theta)
#f = sin(theta)**2

#integral = Integral(f, (theta, 0, pi), (phi, 0, 2*pi)).subs(t, 0)
integral = Integral(f, (phi, 0, 2*pi)).subs(t,0)

