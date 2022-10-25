# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 13:16:28 2022

@author: Pena
"""

from numpy import radians, sin, cos, pi
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d

alphaDeg = 23.5
thetaDeg = 90
phiDeg = 0
wTrans = 2*pi/365
wRot = 732*pi/365
plot = 0
plotArrows = 0

alpha = radians(alphaDeg)
theta = radians(thetaDeg)
phi = radians(phiDeg)
sinAlpha = sin(alpha)
cosAlpha = cos(alpha)


def sun(t=0, alpha=alpha):
    return np.array([
        cos(wRot*t)*cos(wTrans*t) + cosAlpha*sin(wRot*t)*sin(wTrans*t),
        -sin(wRot*t)*cos(wTrans*t) + cosAlpha*cos(wRot*t)*sin(wTrans*t),
        sinAlpha*sin(wTrans*t)])

def place(theta=radians(90), phi=0):
    return np.array([
        sin(theta)*cos(phi),
        sin(theta)*sin(phi),
        cos(theta)+0*phi])

def sunOverPlaceT(t, theta, phi):
    return np.abs(np.inner(sun(t), place(theta, phi)))
def sunOverPlaceTheta(theta, phi, t=0):
    return np.abs(np.inner(sun(t), place(theta, phi)))
def sunOverPlacePhi(phi, theta, t=0):
    return np.abs(np.inner(sun(t), place(theta, phi)))

t = 0
t0 = 0
t1 = 365
phi0 = -pi
phi1 = pi
theta0 = 0
theta1 = pi
#result = integrate.quad(sunOverPlace, t0, t1, args=(theta, phi, S, P),
#                        limit=1000, epsabs=1.49e-3, epsrel=1.49e-3)
#result = integrate.tplquad(sunOverPlace, t0, t1, args=(S, P),
#                        limit=1000, epsabs=1.49e-3, epsrel=1.49e-3)
#result = integrate.quad(sunOverPlacePhi, phi0, phi1, args=(radians(90), 0),
#                        limit=1000, epsabs=1.49e-3, epsrel=1.49e-3)
result = integrate.dblquad(sunOverPlaceTheta, phi0, phi1, theta0, theta1)
#print(result[0]/365)
print(result)

if plot:
   
    ######### Setup PyPlot ##########
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d', proj_type = 'persp')
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)    
    #ax.set_xlim(0, 5000)
    #ax.set_ylim(-2500, 2500)
    #ax.set_zlim(0, 5000)
    #ax.axis("off")
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    ax.set_box_aspect((1,1,1))
    
    
    ######### Plot Sphere #######
    r = 1
    u, v = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
    x = r*sin(u)*cos(v)
    y = r*sin(u)*sin(v)
    z = r*cos(u)
    ax.plot_surface(x, y, z,  rstride=1, cstride=1, color='grey', 
                    alpha=0.1, zorder=2)
    r = 0.5
    ax.plot_surface(r*x, r*y, r*z,  rstride=1, cstride=1, color='grey', 
                    alpha=0.1, zorder=1)
    
    ######### Plot Place on Earth #########
    phi = np.arange(phi0, phi1, 0.1)
    r = 0.5
    px, py, pz = place(radians(80), phi)
    ax.scatter (r*px, r*py, r*pz, color="blue", s=20, zorder= 3)
    if plotArrows:
        ax.quiver(r*px, r*py, r*pz, # <-- starting point of vector
            r*px, r*py, r*pz, # <-- directions of vector
            color = 'blue', alpha = 0.5, lw = 1)
        
    ######### Plot Sun position on Time ########
    t = np.arange(t0, t1, 0.25)
    t = 0
    sx, sy, sz = sun(t)
    ax.scatter(sx, sy, sz, color='red', s=20, zorder=5)
    if plotArrows:
        ax.quiver(r*px, r*py, r*pz, # <-- starting point of vector
            r*sx, r*sy, r*sz, # <-- directions of vector
            color = 'red', alpha = 0.5, lw = 1)
    
    
    plt.show()