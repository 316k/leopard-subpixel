#!/usr/bin/python
# -*- coding: utf-8 -*-
# Visualize checkerboard corners

from __future__ import print_function

import numpy as np
import cv2

from sys import argv

W, H = 9, 6 # Number of x/y intersections to find

# termination criteria
criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)

# prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
objp = np.zeros((H*W,3), np.float32)
objp[:,:2] = np.mgrid[0:W,0:H].T.reshape(-1,2)

for fname in argv[1:]:
    img = cv2.imread(fname)
    gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)

    # Find the chess board corners
    ret, corners = cv2.findChessboardCorners(gray, (W,H), None)

    # If found, add object points, image points (after refining them)
    if ret == True:
        # Draw and display the corners
        cv2.cornerSubPix(gray, corners, (11,11), (-1,-1), criteria)
        cv2.drawChessboardCorners(img, (W,H), corners, ret)

        cv2.imshow(fname, img)

        key = cv2.waitKey(0)
        if key == 27:
            break

        cv2.destroyWindow(fname)

    else:
        print("No match for :", fname)
