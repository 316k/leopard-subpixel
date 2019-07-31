#!/usr/bin/python3

from img import *

import argparse
arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("luts", help="Input phases", nargs='+')
arg_parser.add_argument("--out", help="Merged lut", default='merged.png')

args = arg_parser.parse_args()


imgs = list(map(lambda x: read_img(x)[0], args.luts))

lut = imgs[0].copy()

for test in imgs[1:]:
    assert test.shape == imgs[0].shape, "All luts should have the same shape"

    for i in range(imgs[0].shape[0]):
        for j in range(imgs[0].shape[1]):
            if test[i, j, 2] < lut[i, j, 2]:
                lut[i, j, [0,1]] = test[i, j, [0,1]]


write_img(args.out, lut)

