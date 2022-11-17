# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 14:47:55 2022

@author: Pena
"""


from numpy import radians, sin, cos, pi
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt


r = 1
alphaDeg = 23.5
w = 2*pi/365
W = 732*pi/365
plot = 0
plotArrows = 0

alpha = radians(alphaDeg)
sinAlpha = sin(alpha)
cosAlpha = cos(alpha)

def areaElement(theta, phi, r=1):
    return r**2*sin(theta)

def sun(t=0, alpha=alpha):
    return np.array([
        cos(W*t)*cos(w*t) + cosAlpha*sin(W*t)*sin(w*t),
        -sin(W*t)*cos(w*t) + cosAlpha*cos(W*t)*sin(w*t),
        sinAlpha*sin(w*t)])

def place(theta=radians(90), phi=0):
    return np.array([
        sin(theta)*cos(phi),
        sin(theta)*sin(phi),
        cos(theta)+0*phi])


def sunOverPlace(theta, phi, t):
    return abs((cos(W*t)*cos(w*t) + cosAlpha*sin(W*t)*sin(w*t)) * sin(theta)*cos(phi) + 
    (-sin(W*t)*cos(w*t) + cosAlpha*cos(W*t)*sin(w*t)) * sin(theta)*sin(phi) + 
    sinAlpha*sin(w*t) * cos(theta))*sin(theta)

def sunOverPlaceT(t, theta, phi):
    return np.abs(np.inner(sun(t), place(theta, phi)))*sin(theta)

def sunOverPlaceTheta(theta, phi, t=0):
    return sunOverPlaceT(t, theta, phi)

#def sunOverPlace(t, theta, phi):
#    return np.abs(np.inner(sun(t), place(theta, phi)))*sin(theta)

t0 = 0
t1 = 365
phi0Deg = 0
phi1Deg = 360
theta0Deg = 0
theta1Deg = 180
thetaSteps = 100 * (theta1Deg - theta0Deg)/360 #steps for turn
phiSteps = 100 * (phi1Deg - phi0Deg)/360 #steps for turn
tSteps = 100 * (t1-t0) #steps for day

phi0 = radians(phi0Deg)
phi1 = radians(phi1Deg)
theta0 = radians(theta0Deg)
theta1 = radians(theta1Deg)

theta = pi/2
phi = 0
t = 0
#result = integrate.quad(sunOverPlaceT, t0, t1, args=(theta, phi),
#                        limit=1000, epsabs=1.49e-3, epsrel=1.49e-3)
#result = integrate.quad(sunOverPlaceTheta, theta0, theta1, args=(phi, t),
#                       limit=1000, epsabs=1.49e-3, epsrel=1.49e-3)
result = integrate.dblquad(sunOverPlaceTheta, phi0, phi1, theta0, theta1)
print(result)

phi0 = radians(phi0Deg)
phi1 = radians(phi1Deg)
theta0 = radians(theta0Deg)
theta1 = radians(theta1Deg)

t = np.linspace(t0, t1, int(tSteps))
theta = np.linspace(theta0, theta1, int(thetaSteps))
phi = np.linspace(phi0, phi1, int(phiSteps))
T, THETA, PHI = np.meshgrid(t, theta, phi)

f = sunOverPlace(THETA, PHI, T)
#result = np.trapz(np.trapz(f, phi, axis=0), theta, axis=0)
result = np.trapz(np.trapz(np.trapz(f, t, axis=1), theta, axis=0), phi, axis=0)


print(result)