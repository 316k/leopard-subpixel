#!/usr/bin/python3
# Finds the rigid transform that converts transformed points to
# reference points. Useful for finding the rotation/translation
# between two projectors.

import numpy as np

from sys import argv


def pose_matrix(rotation_mat, translation_mat):
    return np.identity(4).dot(translation_mat.dot(rotation_mat))

def read_transform(fname):
    """Lit les paramètres de la caméra depuis un fichier"""
    f = open(fname, 'r').readlines()

    # Split into lines (lines starting with # are ignored)
    lines = [line.strip() for line in f if len(line.strip()) and line[0] != "#"]

    assert len(lines) == 6, "Unrecognized transformation file: " + fname

    rotation = np.array(list(map(lambda x: list(map(float, x.split('\t'))), lines[0:3])))
    translation = np.array(list(map(float, lines[3:6])))

    translation_mat = np.identity(4)
    translation_mat[:-1, -1] = -translation

    rotation_mat = np.identity(4)
    rotation_mat[:-1, :-1] = rotation

    return pose_matrix(rotation_mat, translation_mat)

if len(argv) != 4:
    print('usage: apply-transform.py transform.dat pts.txt out.txt')
    exit(1)

transform = read_transform(argv[1])

data = np.loadtxt(argv[2])

# Euclidian -> Projective (Nx4)
data = np.hstack((data, np.ones((data.shape[0],1))))

# 4x4 . 4xN -> 4xN .T -> Nx4
data_transformed = transform.dot(data.T).T

# Projective -> Euclidian (Nx3)
data_transformed = data_transformed[:, 0:3] / (data_transformed[:, 3]).reshape((data.shape[0], 1))


np.savetxt(argv[3], data_transformed)
