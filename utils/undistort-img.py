#!/usr/bin/python
# -*- coding: utf-8 -*-
# Undistort an image using a camera's calibration infos

from __future__ import print_function, division

import cv2
import numpy as np
from sys import argv

from helpers import *

if len(argv) < 3:
    print('usage: undistort-img calib.dat img [out-img]')
    exit(1)

mtx, pose, dist = calib_params(argv[1])

img = cv2.imread(argv[2])

out_img = False

if len(argv) == 4:
    out_img = cv2.imread(argv[3])

h, w = img.shape[:2]
newcameramtx, roi = cv2.getOptimalNewCameraMatrix(mtx, dist, (w,h), 1, (w,h))

# undistort
mapx, mapy = cv2.initUndistortRectifyMap(mtx, dist, None, newcameramtx, (w,h), 5)
dst = cv2.remap(img, mapx, mapy, cv2.INTER_LINEAR)

# crop the image
# x,y,w,h = roi
# dst = dst[y:y+h, x:x+w]

if out_img:
    cv2.imwrite(out_img, dst)
else:
    cv2.imshow('Undistorted', dst)

    while True:
        key = cv2.waitKey(10) & 0xFF

        if key == 27:
            break

