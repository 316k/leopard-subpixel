#!/usr/bin/python3
# Finds the rigid transform that converts transformed points to
# reference points. Useful for finding the rotation/translation
# between two projectors.
#
# Outliers are eliminated using a RANDSAC method

import numpy as np
from scipy.linalg import svd
from math import sin, cos

from sys import argv


keep_best = 0.5 # Keep only top 50% points -- score computed with RANDSAC

if len(argv) < 3:
    print('usage: find-rigid-transform.py ref-data.txt transformed-data.txt [keep_best={}]'.format(keep_best))
    exit(1)

ref_data = np.loadtxt(argv[1])
trans_data = np.loadtxt(argv[2])

if len(argv) >= 4:
    keep_best = float(argv[3])

    assert keep_best < 1 and int(keep_best * ref_data.shape[0]) > 0

assert ref_data.shape[0] == trans_data.shape[0]

def normalized(v):
    return v / np.linalg.norm(v)

# https://stackoverflow.com/a/13849249
def vector_angle(v1, v2, opt=False):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    a = normalized(v1)
    b = normalized(v2)

    angle = np.arccos(np.clip(np.dot(a, b), -1.0, 1.0))

    if opt and angle < 0.5: # ~30 deg
        angle = np.linalg.norm(np.arcsin(np.cross(a, b)))

    return angle

def rotation_matrix(angle, rotation_axis):

    rotation_axis = normalized(rotation_axis)
    u, v, w = rotation_axis

    u2, v2, w2 = u*u, v*v, w*w

    cos_a, sin_a = cos(angle), sin(angle)
    icos_a = 1 - cos(angle)

    m1 = u2 + (v2 + w2) * cos_a
    m2 = u * v * icos_a - w * sin_a
    m3 = u * w * icos_a + v * sin_a

    m4 = u * v * icos_a + w * sin_a
    m5 = v2 + (u2 + w2) * cos_a
    m6 = v * w * icos_a - u * sin_a

    m7 = u * w * icos_a - v * sin_a
    m8 = v * w * icos_a + u * sin_a
    m9 = w2 + (u2 + v2) * cos_a

    return np.array([
        [m1, m2, m3],
        [m4, m5, m6],
        [m7, m8, m9],
    ])

def triangle_normal(vertex):
    a, b, c = vertex

    return normalized(np.cross(b - a, c - b));


def transform_dataset(rotation, translation, data):
    return rotation.dot(data.T).T + translation


# Compute scores for the

errors = np.zeros((ref_data.shape[0],))

for _ in range(1000):

    # Pick 3 points
    a = b = c = 0

    while a == b or b == c or c == a:
        a, b, c = np.random.randint(ref_data.shape[0], size=(3,))

    triplet = [a,b,c]

    ref_points = ref_data[triplet]
    trans_points = trans_data[triplet]

    # Compute rotation/translation for those 3 points
    ref_normal = triangle_normal(ref_points)
    trans_normal = triangle_normal(trans_points)

    rotation = rotation_matrix(vector_angle(ref_normal, trans_normal), np.cross(ref_normal, trans_normal))

    translation = (trans_points - rotation.dot(ref_points.T).T).mean(axis=0)

    test_trans_data = transform_dataset(rotation, translation, ref_data)

    # Accumulate error for each point
    errors += np.linalg.norm(test_trans_data - trans_data, axis=1)


# Keep only the best {keep_best} points
# Other points are outliers
threshold = sorted(errors)[int(keep_best * errors.shape[0])]

idx = errors < threshold

ref_data = ref_data[idx]
trans_data = trans_data[idx]


# Compute the least squares rotation/translation for the non-outliers
mean_ref = ref_data.mean(axis=0)

mean_trans = trans_data.mean(axis=0)

ref_centered = ref_data - mean_ref
trans_centered = trans_data - mean_trans

cross_covariance = ref_centered.T.dot(trans_centered)

u, w, vt = svd(cross_covariance, lapack_driver='gesvd')

rotation = u.dot(vt)

rotated_trans = rotation.dot(trans_data.T).T

translation = rotated_trans.mean(axis=0) - mean_ref


def print_arr(arr):
    for line in arr:
        print('\t'.join(map(str, line)))
    print()


print('# Rotation')
print_arr(rotation)

print('# Translation')
print_arr(translation.reshape((3,1)))

recovered_ref_data = (rotated_trans - translation)
print('# Transformation error:', abs(ref_data - recovered_ref_data).sum())
