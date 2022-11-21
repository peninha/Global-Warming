# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 14:47:55 2022

@author: Pena
"""

from numpy import radians, sin, cos, pi
import numpy as np

r = 1
Is = 1
alphaDeg = 23.5
w = 2*pi/365
W = 732*pi/365

alpha = radians(alphaDeg)
sinAlpha = sin(alpha)
cosAlpha = cos(alpha)

def areaElement(theta, r=1):
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

def sunOverPlace(theta, phi, t, r=1):
    return abs((cos(W*t)*cos(w*t) + cosAlpha*sin(W*t)*sin(w*t)) * sin(theta)*cos(phi) + 
    (-sin(W*t)*cos(w*t) + cosAlpha*cos(W*t)*sin(w*t)) * sin(theta)*sin(phi) + 
    sinAlpha*sin(w*t) * cos(theta))   *   r**2*sin(theta)

#############################

"""
r = 1
AreaPx = np.empty(pixels)
for n in range(pixels):
    theta0 = n*delta
    theta1 = (n+1)*delta
    AreaPx[n] = r**2 * delta * ( cos(n*delta) - cos((n+1)*delta) )
np.save("data/AreaPx", AreaPx)
"""
#Area of a each pixel in steradians [units of r^2]
AreaPx = np.load("data/AreaPx.npy")

#############################

"""
r = 1
AreaSlice = np.empty(pixels)
for n in range(pixels):
    theta0 = n*delta
    theta1 = (n+1)*delta
    AreaSlice[n] = r**2 * 2*pi * ( cos(n*delta) - cos((n+1)*delta) )
np.save("data/AreaSlice", AreaSlice)
"""
#Area of a each theta slide in steradians [units of r^2]
AreaSlice = np.load("data/AreaSlice.npy")

#############################

"""
pixels = 512
t0 = 0
t1 = 365
phi0 = 0
phi1 = 2*pi
thetaSteps = 5 #steps from theta0 to theta1
phiSteps = 200 #steps from phi0 to phi1
tSteps = 200 * (t1 - t0) #steps in a 360 degree turn
Is = 1
r = 1

delta = pi/pixels #radians

RadPx = np.empty(pixels)
for n in range(pixels):
    print(n)
    theta0 = n*delta
    theta1 = (n+1)*delta
    t = np.linspace(t0, t1, int(tSteps))
    theta = np.linspace(theta0, theta1, int(thetaSteps))
    phi = np.linspace(phi0, phi1, int(phiSteps))
    T, THETA, PHI = np.meshgrid(t, theta, phi)
    f = sunOverPlace(THETA, PHI, T)
    result = np.trapz(np.trapz(np.trapz(f, t, axis=1), theta, axis=0), phi, axis=0)
    RadPx[n] = Is * result / (2 * (t1-t0)) / AreaSlice[n]
np.save("data/RadPx", RadPx)
"""
#Average Received Radiance [W/m2] over each pixel
RadPx = Is * np.load("data/RadPx.npy")

#############################

#Albedo [0 to 1] from each Pixel 
AlbedoPx = np.load("data/AlbedoPx.npy") 

#Reflected Radiance from each Pixel [W/m2]
ReflexRadPx = AlbedoPx * (RadPx)[:, np.newaxis]

#Average Received Radiance [W/m2] over each theta slice (same as each pixel)
RadSlice = RadPx

#Average Albedo [0 to 1] from each theta slice (average over pixel albedo)
AlbedoSlice = np.sum(AlbedoPx, axis=1) / len(AlbedoPx[0])

#Average Reflected Radiance [W/m2] from each theta slice
ReflexRadSlice = AlbedoSlice * RadSlice

#Average Albedo [0 to 1] of the planet surface
AlbedoSurface = np.sum(ReflexRadSlice*AreaSlice)/np.sum(RadSlice*AreaSlice)

#Fraction [%] of the whole planet Radiation [W] that falls in a single pixel
FracRadiationPx = RadPx*AreaPx / np.sum(RadSlice * AreaSlice)

#Fraction [%] of the whole planet Radiation [W] that falls in a single theta slide
FracRadiationSlice = RadSlice*AreaSlice / np.sum(RadSlice * AreaSlice)

r = 6371 #Earth Radius in km
AreaPxKm = AreaPx * r**2 #Area of each pixel in km2
AreaSliceKm = AreaSlice * r**2 #Area of each theta slide in km2

#Fraction [%] of the whole planet Radiation [W] that falls in a km2
FracRadiationKm = RadSlice / np.sum(RadSlice * AreaSliceKm)

AlbedoAtmos = 0.232
AbsorbedAtmos = 0.298
AlbedoPlanet = AlbedoSurface*(1-AlbedoAtmos)*(1-AbsorbedAtmos) + AlbedoAtmos
