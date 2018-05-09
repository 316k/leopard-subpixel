#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division, print_function

from sys import argv, stdout
import numpy as np

import cv2

from img import *

img, depth = read_img(argv[1])
out = np.zeros((img.shape[0], img.shape[1]))

i=0
for y in range(img.shape[0]):
    for x in range(img.shape[1]):
        out[y][x] = 255 if img[y][x][2]/65535 < 0.2 else 0
        if img[y][x][2]/65535 < 0.2:
            i += 1

print(i)
write_img("mask.pgm", out, 8)

