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

ColorPx = np.empty(shape=(height, width), dtype='<U6')
AlbedoPx = np.zeros(shape=(height, width))
for n in range(height):      # this row
    for m in range(width):   # and this row was exchanged
        r, g, b = pixels[m, n]
        ColorPx[n, m] = (f"{r:02x}{g:02x}{b:02x}")
        AlbedoPx[n, m] = Color2Albedo[ColorPx[n, m]]

np.save("data/ColorPx", ColorPx)
np.save("data/AlbedoPx", AlbedoPx)