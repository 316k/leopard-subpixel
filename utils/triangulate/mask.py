#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filter pixels according to the Blue channel

from __future__ import division, print_function

from sys import argv, stdout
import numpy as np

import cv2

from img import *

threshold = 0.2

if len(argv) < 3:
    print('usage: mask.py img.ppm out.ppm [threshold=%f]' % threshold)
    exit(1)

img, depth = read_img(argv[1])

out_fname = argv[2]

if len(argv) >= 4:
    threshold = float(argv[3])

out = np.zeros(img.shape)

i=0
for y in range(img.shape[0]):
    for x in range(img.shape[1]):

        if img[y][x][2]/65535 < threshold:
            out[y][x] = img[y][x]
            i += 1
        else:
            out[y][x][0] = out[y][x][1] = out[y][x][2] = 65535

print(i, '/', img.shape[0] * img.shape[1], '(', i / (img.shape[0] * img.shape[1]), ')')

write_img(out_fname, out, 16)

