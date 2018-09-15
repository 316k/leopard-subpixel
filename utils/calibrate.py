#!/usr/bin/python
# -*- coding: utf-8 -*-
# Calibrate a camera from checkerboards images
# Code adapted from the OpenCV documentation

from __future__ import print_function, division

import numpy as np
import cv2

from sys import argv

from helpers import eprint

# Args
arg0 = argv[0]

argv = argv[1:]

verbose = False
no_radial_disto = False
external_from_first_checkerboard = False

while len(argv) and argv[0].startswith('-') and argv[0] != '--':
    if argv[0] == '-v':
        verbose = True
    elif argv[0] == '-p':
        no_radial_disto = True
    elif argv[0] == '-e':
        external_from_first_checkerboard = True

    argv = argv[1:]

if len(argv) == 0:
    print('usage: calibrate.py [-v verbose] [-p no radial disto] [-e] checkerboards...')
    print('\t-e: pick rotation/translation from the first checkerboard')
    exit(1)

# Actual stuff
W, H = 9, 6

# termination criteria
criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)

# prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
objp = np.zeros((H*W,3), np.float32)
objp[:,:2] = np.mgrid[0:W,0:H].T.reshape(-1,2)

# Arrays to store object points and image points from all the images.
objpoints = [] # 3d point in real world space
imgpoints = [] # 2d points in image plane.

valid_images = []

for i, fname in enumerate(argv):
    img = cv2.imread(fname)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    # Find the chess board corners
    ret, corners = cv2.findChessboardCorners(gray, (W,H), None)

    # If found, add object points, image points (after refining them)
    if verbose:
        print(fname, ":", ret)

    if ret == True:
        objpoints.append(objp)

        cv2.cornerSubPix(gray,corners,(11,11),(-1,-1),criteria)

        imgpoints.append(corners)

        valid_images.append(fname)
    elif external_from_first_checkerboard and i == 0:
        eprint('Error: First checkerboard not detected')
        exit(1)

flags = []

# Force la distortion radiale à zéro
if no_radial_disto:
    flags=[
        cv2.CALIB_FIX_K1, cv2.CALIB_FIX_K2, cv2.CALIB_FIX_K3,
        cv2.CALIB_FIX_K4, cv2.CALIB_FIX_K5, cv2.CALIB_FIX_K6
    ]

ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(objpoints, imgpoints, gray.shape[::-1], None, None, flags=reduce(lambda acc, x: acc | x, flags, 0))

def print_arr(arr):
    for line in arr:
        print('\t'.join(map(str, line)))
    print()

print("# Internes")
print_arr(mtx)

print("# Rotations")
for i, r in enumerate(rvecs):
    print('#', valid_images[i])
    print_arr(cv2.Rodrigues(r)[0])

    if external_from_first_checkerboard:
        break

print("# Translations")
for i, t in enumerate(tvecs):
    print('#', valid_images[i])
    print_arr(t)

    if external_from_first_checkerboard:
        break

print("# Distortion coeffs")
print_arr(dist)

mean_error = 0
for i in xrange(len(objpoints)):
    imgpoints2, _ = cv2.projectPoints(objpoints[i], rvecs[i], tvecs[i], mtx, dist)
    error = cv2.norm(imgpoints[i],imgpoints2, cv2.NORM_L2)/len(imgpoints2)
    mean_error += error

if verbose:
    print("total error: ", mean_error/len(objpoints))
