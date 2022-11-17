# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 00:49:56 2022

@author: Pena
"""

import numpy as np
import pickle
from PIL import Image
from numpy import sin, cos, pi
with open('data/Color2Albedo.pkl', 'rb') as f:
    Color2Albedo = pickle.load(f)


img = Image.open('imgs/SurfaceAlbedo_2021.png')
pixels = img.load() 
width, height = img.size


Hex = np.empty(shape=(height, width), dtype='<U6')
Albedo = np.zeros(shape=(height, width))
AreaWeight = np.zeros(height)
for x in range(height):      # this row
    for y in range(width):   # and this row was exchanged
        r, g, b = pixels[y, x]
        Hex[x, y] = (f"{r:02x}{g:02x}{b:02x}")
        #Albedo[x, y] = Color2Albedo[Hex[x,y]]
    AreaWeight[x] = pi/512*(cos(x*pi/512) - cos((x+1)*pi/512))
        # in case your image has an alpha channel
        # r, g, b, a = pixels[x, y]

        #print(x, y, f"#{r:02x}{g:02x}{b:02x}")
   
