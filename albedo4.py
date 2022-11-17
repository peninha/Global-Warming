# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 14:47:55 2022

@author: Pena
"""

from numpy import radians, degrees, sin, cos, pi
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

r = 1
alphaDeg = 23.5
w = 2*pi/365
W = 732*pi/365

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


t0 = 0
t1 = 365
phi0Deg = -180
phi1Deg = 180
#theta0Deg = 80
#theta1Deg = 100
pixels = 512
thetaSteps = 11 #steps from theta0 to theta1
steps = 50 #steps in a 360 degree turn

phi0 = radians(phi0Deg)
phi1 = radians(phi1Deg)
#theta0 = radians(theta0Deg)
#theta1 = radians(theta1Deg)
delta = pi/pixels #radians
deltaDegree = degrees(delta)

"""
FluxoMedio = np.empty(pixels)
for i in range(148, 169):
    print(i)
    theta0 = i*delta
    theta1 = (i+1)*delta
    phiSteps = steps * (phi1Deg - phi0Deg)/360 #steps for turn
    tSteps = steps * (t1-t0) #steps for day
    t = np.linspace(t0, t1, int(tSteps))
    theta = np.linspace(theta0, theta1, int(thetaSteps))
    phi = np.linspace(phi0, phi1, int(phiSteps))
    T, THETA, PHI = np.meshgrid(t, theta, phi)
    f = sunOverPlace(THETA, PHI, T)
    result = np.trapz(np.trapz(np.trapz(f, t, axis=1), theta, axis=0), phi, axis=0)/2
    FluxoMedio[i] = result / (t1-t0)

"""
FluxoMedio = np.load("data/FluxoMedioT11S50.npy")

Area = np.empty(pixels)
for i in range(pixels):
    theta0 = i*delta
    theta1 = (i+1)*delta
    superficie = integrate.dblquad(areaElement, phi0, phi1, theta0, theta1)
    Area[i] = superficie[0]

DensidadeFluxoMedio = FluxoMedio/Area
totalArea = np.sum(Area)
totalFluxoMedio = np.sum(FluxoMedio)
totalDensidadeFluxoMedio = totalFluxoMedio/totalArea
