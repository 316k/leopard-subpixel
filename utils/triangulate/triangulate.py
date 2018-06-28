#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division, print_function

from sys import argv, stdout, stderr
import numpy as np
import cv2
from colorsys import hsv_to_rgb

from img import *

def pose_matrix(rotation_mat, translation_mat):
    return np.identity(4)[:-1, :].dot(rotation_mat.dot(translation_mat))

def params(fname):
    """Lit les paramètres de la caméra depuis un fichier"""
    f = open(fname, 'r').readlines()

    # Découpe en lignes (lignes qui commencent par "#" == commentaire)
    lines = [line.strip() for line in f if len(line.strip()) and line[0] != "#"]

    assert len(lines) == 10

    internes = np.array(map(lambda x: map(float, x.split('\t')), lines[0:3]))
    rotation = np.array(map(lambda x: map(float, x.split('\t')), lines[3:6]))
    translation = np.array(map(float, lines[6:9]))
    distortion_coeffs = np.array(map(float, lines[9].split('\t')))

    translation_mat = np.identity(4)
    translation_mat[:-1, -1] = translation

    rotation_mat = np.identity(4)
    rotation_mat[:-1, :-1] = rotation

    return internes, pose_matrix(rotation_mat, translation_mat), distortion_coeffs

# triangulation::matrixCorr
def compute_pixel_points(lut1, lut2, threshold=0.1):

    ref_h, ref_w = lut1.shape[0], lut1.shape[1]
    matching_h, matching_w = lut2.shape[0], lut2.shape[1]

    good_points = lut1[lut1[:, :, 2]/65535 < threshold]

    ref_pts = np.zeros((2, good_points.shape[0]))
    matching_pts = np.zeros((2, good_points.shape[0]))

    i = 0
    for y, row in enumerate(lut1): # TODO fold en un seul np.ndenumerate
        for x, col in enumerate(row):
            # col = [x, y, error]
            if col[2]/65535 < threshold:
                ref_pts[0][i] = x
                ref_pts[1][i] = y
            
                matching_pts[0][i] = col[0]/65535 * matching_w
                matching_pts[1][i] = col[1]/65535 * matching_h

                i += 1

    assert i == good_points.shape[0]

    return ref_pts, matching_pts

# undistortMatrix
def undistort_matrix(points, internes, disto):

    n = points.shape[1]

    points1D = points.transpose().reshape(points.shape[1], 1, 2)

    undistorted_points = cv2.undistortPoints(points1D, internes, disto);

    return undistorted_points.reshape((n, 2)).transpose()

# lut2corr
def lut_to_corrected_points(lut1, internes1, disto1, lut2, internes2, disto2, threshold=0.1):

    ref_pts, matching_pts = compute_pixel_points(lut1, lut2, threshold)

    undistorted_points1 = undistort_matrix(ref_pts, internes1, disto1)
    undistorted_points2 = undistort_matrix(matching_pts, internes2, disto2)

    return undistorted_points1, undistorted_points2


if len(argv) < 5:
    print("Usage %s lut1.ppm lut2.ppm cam1params.dat cam2params.dat [lut3d.ppm] > pts.txt" % argv[0])
    print()
    print("  Outputs 3D points from two LUT + calibration data on stdout")
    print("  If a last ppm filename is given, writes a LUT with RGB -> (x, y, z)")
    exit(-1)

# Load data
internes1, pose1, disto1 = params(argv[3])
internes2, pose2, disto2 = params(argv[4])

lut1, depth = read_img(argv[1])
assert depth == 65535

lut2, depth = read_img(argv[2])
assert depth == 65535

lut3d = False
threshold = 0.40

# Dump une depth map en ppm si demandé
if len(argv) > 5:
    lut3d = np.zeros((lut1.shape[0], lut1.shape[1], 3))

    grid = np.zeros((lut1.shape[0], lut1.shape[1], 2), dtype=int)

    for y, row in enumerate(lut1):
        for x, col in enumerate(row):
            grid[y][x][0] = x
            grid[y][x][1] = y

    xy = grid[lut1[:, :, 2]/65535 < threshold]

# Normalize/undistort lut points
undistorted_points1, undistorted_points2 = lut_to_corrected_points(lut1, internes1, disto1, lut2, internes2, disto2, threshold)

# Triangulate
pts = cv2.triangulatePoints(pose1, pose2, undistorted_points1, undistorted_points2) # output une 4xN

assert pts.shape[0] == 4

# projectif -> euclidien (Nx3)
euclidien = (pts[0:3, :] / np.array([pts[3, :], pts[3, :], pts[3, :]])).transpose()
assert euclidien.shape[1] == 3


minv = np.array([-7, -8, 0])
maxv = np.array([20, 8, 22])

for n, row in enumerate(euclidien):
    if all(row > minv) and all(row < maxv):
        print(row[0], row[1], row[2])

        if lut3d is not False:
            x, y = xy[n]
            lut3d[y][x] = (row - minv)/(maxv - minv) * 65535

if lut3d is not False:
    write_img(argv[5], lut3d, 16)

#np.savetxt(stdout, euclidien, fmt='%.10f %.10f %.10f')
