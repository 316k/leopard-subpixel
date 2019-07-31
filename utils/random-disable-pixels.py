#!/usr/bin/python3

from img import *
import numpy as np

import argparse
arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("seed", help="Random mask seed", type=int, default=1)
arg_parser.add_argument("img", help="Input image")
arg_parser.add_argument("out", help="Maksed image")
arg_parser.add_argument("--block-size", help="Size of a cut out block", type=int, default=5)

args = arg_parser.parse_args()

bs=args.block_size

img, depth = read_img(args.img)

np.random.seed(args.seed)

offsets = np.random.randint(0, bs, (2,))

h, w = img.shape

print(w, h)

low_res_mask = np.random.randint(0, 2, (h // bs + 1, w // bs + 1))

full_mask = np.zeros((h, w))

for y in range(h):
    for x in range(w):
        full_mask[y][x] = low_res_mask[(y - offsets[0]) // bs, (x - offsets[1]) // bs]

write_img(args.out, img * full_mask, depth=16 if depth == 65535 else 8)

