#!/usr/bin/python
# -*- coding: utf-8 -*-
# Capture synchronized patterns
# Display & capture patterns from the same computer

from __future__ import print_function, division

import cv2
import numpy as np

from sys import argv
from glob import glob

from time import time

if len(argv) < 2:
    print('usage: capture-slow.py path/to/imgs')
    exit(1)

ref = sorted(glob(argv[1] + 'leo_*.pgm'))

delay = 0.5
last_capture_time = 0

fps = 30
width = 1280
height = 720

# capture from camera at location 0
cap = cv2.VideoCapture(0)

cap.set(3, width) # Width
cap.set(4, height)  # Height
cap.set(5, fps)   # FPS

# Wonky ?
cap.set(21, 0)   # Auto_exposure
cap.set(39, 0)   # Autofocus

# cap.set(cv2.CAP_PROP_AUTOFOCUS, 0)
# cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, 0)

i = 0
# threshold = 0
save_img = False

# def change_thresh(x):
#     global threshold
#     threshold = x

def set_save(x):
    global save_img, last_capture_time
    save_img = bool(x)
    last_capture_time = time()

cv2.namedWindow('cam')
# cv2.createTrackbar('thresh', 'cam', 0, 255, change_thresh)
cv2.createTrackbar('record', 'cam', 0, 1, set_save)

# References
cv2.namedWindow('ref')
cv2.moveWindow('ref', 0, 0)
current_ref = cv2.imread(ref[i])
cv2.imshow('ref', current_ref)

last = None

while True:

    ret, img = cap.read()
    img = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)

    if last is None:
        last = img
        continue

    # diff = ((img - last)**2 / img.size).sum()
    # print(diff)

    current_time = time()

    if save_img and last_capture_time + delay < current_time:
        # and diff > threshold:
        last_capture_time = current_time
        cv2.imwrite("%03d.pgm" % i, img)

        i += 1

        # Fin de la séquence
        if i == len(ref):
            break

        current_ref = cv2.imread(ref[i])
        cv2.imshow('ref', current_ref)

    if not save_img:
        cv2.imshow("cam", img)

    last = img

    key = cv2.waitKey(10)

    if key == 27:
        break

    if key == 102:
        print(cv2.getWindowProperty("ref", cv2.WND_PROP_FULLSCREEN))
        cv2.setWindowProperty("ref", cv2.WND_PROP_FULLSCREEN, 4);

cv2.destroyAllWindows()
cv2.VideoCapture(0).release()
