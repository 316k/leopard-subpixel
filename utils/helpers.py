#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import cv2
import numpy as np

from sys import stderr

def eprint(*args, **kwargs):
    print(*args, file=stderr, **kwargs)

def pose_matrix(rotation_mat, translation_mat):
    return np.identity(4)[:-1, :].dot(translation_mat.dot(rotation_mat))

def calib_params(fname):
    """Lit les paramètres de la caméra depuis un fichier"""
    f = open(fname, 'r').readlines()

    # Découpe en lignes (lignes qui commencent par "#" == commentaire)
    lines = [line.strip() for line in f if len(line.strip()) and line[0] != "#"]

    assert len(lines) == 10, "Unrecognized calibration file: " + fname

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

def triangulate(lut1, internes1, disto1, pose1, lut2, internes2, disto2, pose2, threshold):

    # Normalize/undistort lut points
    undistorted_points1, undistorted_points2 = lut_to_corrected_points(lut1, internes1, disto1, lut2, internes2, disto2, threshold)

    # Triangulate
    pts = cv2.triangulatePoints(pose1, pose2, undistorted_points1, undistorted_points2) # output une 4xN

    assert pts.shape[0] == 4

    # Projective -> Euclidian (Nx3)
    euclidienPts = (pts[0:3, :] / np.array([pts[3, :], pts[3, :], pts[3, :]])).transpose()
    assert euclidienPts.shape[1] == 3

    return euclidienPts
