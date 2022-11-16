# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 00:49:56 2022

@author: Pena
"""

import numpy as np
import pickle
from PIL import Image

img = Image.open('imgs/SurfaceAlbedo_2021_Blues_Scale.png')
pixels = img.load() 
width, height = img.size

hex = list()
for y in range(height):      # this row
    for x in range(width):   # and this row was exchanged
        r, g, b = pixels[x, y]
        hex.append(f"{r:02x}{g:02x}{b:02x}")
        # in case your image has an alpha channel
        # r, g, b, a = pixels[x, y]

        #print(x, y, f"#{r:02x}{g:02x}{b:02x}")
        
HexColorScale = list()
previous, new = (0, 0)
for i in hex:
    if i == previous:
        new = i
    else:
        if previous != 0 and new != 0:
            HexColorScale.append(new)
            new = 0
    previous = i
if new != 0:
    HexColorScale.append(new)

HexColorScale.reverse()
np.save("data/HexColorScale", HexColorScale)

AlbedoScale = [
0.050,
0.062,
0.075,
0.087,
0.100,
0.112,
0.125,
0.137,
0.150,
0.162,
0.175,
0.187,
0.200,
0.213,
0.225,
0.238,
0.250,
0.263,
0.275,
0.288,
0.300,
0.313,
0.325,
0.338,
0.351,
0.363,
0.376,
0.388,
0.401,
0.413,
0.426,
0.438,
0.451,
0.463,
0.476,
0.489,
0.501,
0.514,
0.526,
0.539,
0.551,
0.564,
0.576,
0.589,
0.601,
0.614,
0.626,
0.639,
0.652,
0.664,
0.677,
0.689,
0.702,
0.714,
0.727,
0.739,
0.752,
0.764,
0.777,
0.790,
0.802,
0.815,
0.827,
]
np.save("data/AlbedoScale", AlbedoScale)

Color2Albedo = dict()
for i in range(len(AlbedoScale)):
    Color2Albedo[HexColorScale[i]] = AlbedoScale[i]

with open('data/Color2Albedo.pkl', 'wb') as f:
    pickle.dump(Color2Albedo, f)