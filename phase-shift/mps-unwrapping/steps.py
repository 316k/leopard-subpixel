#!/usr/bin/python3

from math import ceil, log
from sys import argv

argv = argv[1:]

base = float(argv[0])
low = float(argv[1])
high = float(argv[2])

r = round

def round(x):
    """390.96299999999997 => 390.963"""
    return r(x * 10**4) / 10**4

print(' '.join([str(round(low * base ** i)) for i in range(0, ceil(log(high/low, base)) + 1)]))
