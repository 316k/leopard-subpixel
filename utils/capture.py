#!/usr/bin/python
# -*- coding: utf-8 -*-
# Capture patterns & (try to) extract only non-mixed frames

from __future__ import print_function, division

import cv2
import numpy as np

fps = 16
width = 1280
height = 720

# capture from camera at location 0
cap = cv2.VideoCapture(1)

cap.set(3, width) # Width
cap.set(4, height)  # Height
cap.set(5, fps)   # FPS

# Wonky ?
cap.set(21, 0)   # Auto_exposure
cap.set(39, 0)   # Autofocus

# cap.set(cv2.CAP_PROP_AUTOFOCUS, 0)
# cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, 0)

i = 0
threshold = 0
save_img = False

def change_thresh(x):
    global threshold
    threshold = x


def set_save(x):
    global save_img
    save_img = bool(x)
    
cv2.namedWindow('input')
cv2.createTrackbar('thresh', 'input', 0, 255, change_thresh)
cv2.createTrackbar('record', 'input', 0, 1, set_save)


last = None

while True:

    ret, img = cap.read()
    img = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)

    if last is None:
        last = img
        continue

    diff = ((img - last)**2 / img.size).sum()

    print(diff)

    if save_img and diff > threshold:
        cv2.imwrite("%03d.pgm" % i, img)
        i += 1

    cv2.imshow("input", img)

    last = img

    key = cv2.waitKey(10) & 0xFF

    if key == 27:
        break

cv2.destroyAllWindows()
cv2.VideoCapture(0).release()
