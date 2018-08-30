#!/usr/bin/python
# -*- coding: utf-8 -*-
# Preview a camera with OpenCV
# Type s to take a screenshot

from __future__ import print_function, division

import cv2
import numpy as np

from sys import argv

from glob import glob

fps = 30

camera = 1

if len(argv) > 1:
    camera = int(argv[1])

print("Using camera ", argv)

# capture from camera
cap = cv2.VideoCapture(camera)

cap.set(3, 1280) # Width
cap.set(4, 720)  # Height
cap.set(5, fps)   # FPS

# Wonky ?
cap.set(21, 0)   # Auto_exposure
cap.set(39, 0)   # Autofocus

# cap.set(cv2.CAP_PROP_AUTOFOCUS, 0)
# cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, 0)

i = len(glob('screenshot*.pgm'))
threshold = 12

last = None

while True:

    ret, img = cap.read()
    img = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)

    if last is None:
        last = img
        continue

    diff = ((img - last)**2 / img.size).sum()

    print(diff)

    
    overlay = img.copy()
    font = cv2.FONT_HERSHEY_PLAIN
    cv2.putText(overlay, 'screenshots: ' + str(i), (25,25), font, 1, 200 if i % 2 else 80, 1)

    cv2.imshow("input", overlay)
    
    key = cv2.waitKey(10) & 0xFF

    if key == 27:
        break
    elif key == ord('s'):
        cv2.imwrite("screenshot-%03d.pgm" % i, img)
        i += 1

cv2.destroyAllWindows()
cv2.VideoCapture(0).release()
