#!/usr/bin/python3

from img import *
import numpy as np

import argparse
arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("seed", help="Random mask seed", type=int, default=1)
arg_parser.add_argument("img", help="Input image")
arg_parser.add_argument("out", help="Maksed image")
arg_parser.add_argument("--block-size", help="Size of a cut out block", type=int, default=5)
arg_parser.add_argument("--white", help="Mask with white instead of black", action='store_true')

args = arg_parser.parse_args()

bs=args.block_size

img, depth = read_img(args.img)

np.random.seed(args.seed)

offsets = np.random.randint(0, bs, (2,))

h, w = img.shape[:2]

low_res_mask = np.random.randint(0, 2, (h // bs + 1, w // bs + 1))

full_mask = np.zeros(img.shape)

for y in range(h):
    for x in range(w):
        if len(img.shape) == 2:
            full_mask[y, x] = low_res_mask[(y - offsets[0]) // bs, (x - offsets[1]) // bs]
        else:
            full_mask[y, x, :] = low_res_mask[(y - offsets[0]) // bs, (x - offsets[1]) // bs]

out_img = img * full_mask

if args.white:
    out_img[full_mask == 0] = depth

write_img(args.out, out_img, depth=16 if depth == 65535 else 8)

