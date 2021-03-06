#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division, print_function

#from PIL import Image
import numpy as np
import cv2

# Pillow PNG stuff
def img(path):
    im = Image.open(path)
    return np.array(im.getdata()).reshape((im.width, im.height, 3))

def write(fname, arr):
    im = Image.new('L', arr.shape[0:2])

    out = (arr * 2**8).flatten().tolist()

    im.putdata(out)

    im.save(fname)


# ppm/pgm stuff
def read_pgm(f):
    """Return a raster of integers from a PGM as a list of lists."""
    line = f.readline()

    while line[0] in (35, '#'): # python2/3 shenanigans
        line = f.readline()

    width, height = [int(i) for i in line.split()]

    depth = int(f.readline())

    assert depth == 255 or depth == 65535

    raster = np.zeros((height, width))

    for y in range(height):
        for x in range(width):
            if depth == 255:
                raster[y][x] = ord(f.read(1))
            else:
                raster[y][x] = (ord(f.read(1)) << 8) + ord(f.read(1))

    return raster, depth

def read_ppm(f):
    """Return a raster of integers from a PGM as a list of lists."""
    line = f.readline()

    while line[0] in (35, '#'): # python2/3 shenanigans
        line = f.readline()

    width, height = [int(i) for i in line.split()]

    depth = int(f.readline())

    assert depth == 255 or depth == 65535

    raster = np.zeros((height, width, 3))

    for y in range(height):
        for x in range(width):
            for c in range(3):
                if depth == 255:
                    raster[y][x][c] = ord(f.read(1))
                else:
                    raster[y][x][c] = (ord(f.read(1)) << 8) + ord(f.read(1))

    return raster, depth

def read_img(fname):

    img = cv2.imread(fname, -1)

    if len(img.shape) == 3:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

    return img.astype('float64'), 65535 if img.dtype == 'uint16' else 255

    # f = open(fname, 'rb')
    # header = f.readline()

    # if header == b'P5\n':
    #     return read_pgm(f)
    # elif header == b'P6\n':
    #     return read_ppm(f)

    # print(fname + ": not a ppm or pgm")
    # exit(-1)

def save_pgm(fname, out, depth):

    # dtype = 'u1' if depth == 8 else '>u2'
    # copy = np.array(out, dtype=dtype)
    # cv2.imwrite(fname, copy)

    h, w = out.shape[0], out.shape[1]

    with open(fname, 'wb') as f:

        f.write("P5\n{} {}\n{}\n".format(w, h, (1<<depth) - 1))

        dtype = 'u1' if depth == 8 else '>u2'

        copy = np.array(out, dtype=dtype)

        f.write(copy.tostring())

# TODO : deprecate
def save_ppm(fname, out, depth):

    # dtype = 'u1' if depth == 8 else '>u2'
    # copy = np.array(out, dtype=dtype)
    # #cv2.imwrite(fname, cv2.cvtColor(copy, cv2.COLOR_RGB2BGR))
    # cv2.imwrite(fname, copy)

    h, w = out.shape[0], out.shape[1]

    with open(fname, 'wb') as f:

        f.write("P6\n{} {}\n{}\n".format(w, h, (1<<depth) - 1))

        dtype = 'u1' if depth == 8 else '>u2'

        copy = np.array(out, dtype=dtype)

        f.write(copy.tostring())

def write_img(fname, arr, depth=16):


    img = arr.astype('uint' + str(depth))

    if len(img.shape) == 3:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

    cv2.imwrite(fname, img)

    # if len(arr.shape) == 3:
    #     save_ppm(fname, arr, depth)
    # elif len(arr.shape) == 2:
    #     save_pgm(fname, arr, depth)
    # else:
    #     raise ValueError("Wrong array")

if __name__ == "__main__":
    img = np.zeros((600, 100))

    for i in range(600):
        img[i, 0:50] = i * 100

    write_img("waaaat.pgm", img, 16)
