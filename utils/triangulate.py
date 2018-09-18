#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division, print_function

from sys import argv, stdout
import numpy as np
import cv2
from colorsys import hsv_to_rgb

from img import *
from helpers import *


if len(argv) < 5:
    print("Usage triangulate.py lutProj.ppm lutCam.ppm \\")
    print("\tprojCalib.dat camCalib.dat [lut3d.ppm] > pts.txt")
    print()
    print("  Outputs 3D points from two LUT + calibration data on stdout")
    print("  If a last ppm filename is given, writes a LUT with RGB -> (x, y, z)")
    exit(-1)

# Load data
internes1, pose1, disto1 = calib_params(argv[3])
internes2, pose2, disto2 = calib_params(argv[4])

lut1, depth = read_img(argv[1])
assert depth == 65535

lut2, depth = read_img(argv[2])
assert depth == 65535

lut3d = False
threshold = 1 #0.40

# Export a ppm xyz-map if required
if len(argv) > 5:
    lut3d = np.zeros((lut1.shape[0], lut1.shape[1], 3))

    xy = generate_xy_list(lut1, threshold)

pts = triangulate(lut1, internes1, disto1, pose1, lut2, internes2, disto2, pose2, threshold)

# -- Generate pts + lut3D --
mean = pts.mean(axis=0)
stddev = pts.std(axis=0)

minv = mean - 3 * stddev
maxv = mean + 3 * stddev

eprint("Mean:", mean)
eprint("Stddev:", stddev)
eprint("Min:", minv)
eprint("Max:", maxv)

# minv = np.array([-1e10, -1e10, -1e10])
# maxv = np.array([1e10, 1e10, 1e10])

for n, point in enumerate(pts):
    if all(point > minv) and all(point < maxv):
        print(point[0], point[1], point[2])

        if lut3d is not False:
            x, y = xy[n]
            lut3d[y][x] = (point - minv)/(maxv - minv) * 65535

if lut3d is not False:
    write_img(argv[5], lut3d, 16)

# np.savetxt(stdout, pts, fmt='%.10f %.10f %.10f')
