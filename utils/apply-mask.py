#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filter pixels according to the Blue channel

from __future__ import division, print_function

from sys import argv, stdout
import numpy as np

import cv2

from img import *

def eprint(*args, **kwargs):
    print(*args, file=stderr, **kwargs)

if len(argv) < 4:
    print('usage: apply-mask.py mask.pgm img.ppm out.ppm')
    exit(1)

mask, _ = read_img(argv[1])
img, depth = read_img(argv[2])

assert img.shape[0:2] == mask.shape[0:2], "Image and mask should have the same dimensions"

out_fname = argv[3]

out = img.copy()

for y in range(img.shape[0]):
    for x in range(img.shape[1]):

        if mask[y][x] == 0:
            out[y][x][0] = out[y][x][1] = out[y][x][2] = 65535

write_img(out_fname, out, 16)

