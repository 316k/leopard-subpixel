#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import cv2
import numpy as np

def pose_matrix(rotation_mat, translation_mat):
    return np.identity(4)[:-1, :].dot(rotation_mat.dot(translation_mat))

def calib_params(fname):
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
