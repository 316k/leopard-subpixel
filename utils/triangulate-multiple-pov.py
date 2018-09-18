#!/usr/bin/python3
# Triangulate points from multiple camera point of views
# The projector should stay at the same position during each scan and
# the camera should be the same (or have the same internal parameters
# matrix) for every pair of LUTs

# FIXME : The coordinate (0, 0, 0) is used to specify "no point" in
# lut3d, another "special" value should be used since (0,0,0) is a
# perfectly valid point (something like NaN or INF might do it).


import numpy as np

from sys import argv, stdout

from img import *
from helpers import *

debug = False

if len(argv) > 1 and argv[1] == '-d':
    debug = True
    argv = argv[0] + argv[2:]

if len(argv) < 5:
    print('usage: [-d] triangulate-multiple-pov.py cam.dat proj.dat lutCam1.ppm lutProj1.ppm \\')
    print('\t [lutCam2.ppm lutProj2.ppm [...]] [lut3d.ppm] > pts.txt')
    exit(1)

threshold = 1

cam_internes, cam1_pose, cam_disto = calib_params(argv[1])
proj_internes, proj_pose, proj_disto = calib_params(argv[2])

lut_cam, depth = read_img(argv[3])
assert depth == 65535

lut_proj, depth = read_img(argv[4])
assert depth == 65535


argv = argv[5:]

# Last arg == lut3d (if odd number of args)
if len(argv) % 2 == 1:
    lut3d_fname = argv[-1]
    argv = argv[:-1]
else:
    lut3d_fname = False

assert len(argv) % 2 == 0

cam_luts_fnames = argv[::2]
proj_luts_fnames = argv[1::2]


# Compute the initial lut3d
lut3d = np.zeros(lut_proj.shape)

xy = generate_xy_list(lut_proj, threshold)

pts = triangulate(lut_proj, proj_internes, proj_disto, proj_pose, lut_cam, cam_internes, cam_disto, cam1_pose, threshold)

h, w = lut3d.shape[0], lut3d.shape[1]

mean = pts.mean(axis=0)
stddev = pts.std(axis=0)

minv = mean - 3 * stddev
maxv = mean + 3 * stddev


pts[(pts > maxv).any(axis=1)] = (0, 0, 0)
pts[(pts < minv).any(axis=1)] = (0, 0, 0)

# Blue channel of the LUT == match distance (lower is better)
confidence = np.ones(lut_proj.shape[0:2]) * float('inf')

for n, point in enumerate(pts):
    x, y = xy[n]
    lut3d[y][x] = point

    if not (point == 0).all():
        confidence[y, x] = lut_proj[y, x, 2]

def save_lut3d(i):

    out = ((lut3d.reshape((h * w, 3)) - minv)/(maxv - minv)).reshape((h, w, 3)) * 65535

    out[(lut3d == 0).all(axis=2)] = (0, 0, 0)

    write_img('lut3d-' + str(i) + '.ppm', out, 16)


if debug:
    save_lut3d(0)


# Iteratively upgrade lut3d by using matches between lut3d and current lut
i = 0
for lut_cam_fname, lut_proj_fname in zip(cam_luts_fnames, proj_luts_fnames):

    i += 1

    eprint('#', lut_cam_fname, lut_proj_fname)
    
    lut_cam, depth = read_img(lut_cam_fname)
    assert depth == 65535

    lut_proj, depth = read_img(lut_proj_fname)
    assert depth == 65535

    # Find projector pixels seen before (at any iteration) & now
    lut3d_mask = np.zeros(lut3d.shape[0:2], dtype='uint8')
    lut3d_mask[(lut3d > 0).any(axis=2)] = 1

    lut_proj_mask = np.zeros(lut_proj.shape[0:2], dtype='uint8')
    lut_proj_mask[lut_proj[:, :, 2] < threshold * 65535] = 1

    intersection = lut3d_mask & lut_proj_mask

    # Find new camera pose

    # World points seen before & now : Nx3 points (x,y,z)
    world_pts = lut3d[intersection > 0]

    assert len(world_pts.shape) == 2 and world_pts.shape[1] == 3

    # Array of 2D points (projector pixels matching the world points)
    img_pts_proj = np.argwhere(intersection > 0)

    # Array of 2D points (camera pixels matching the world points)
    img_pts = np.zeros(img_pts_proj.shape, dtype='float64')

    for n, pt in enumerate(img_pts_proj):

        y, x = pt

        img_pts[n][0] = lut_proj[y][x][0] / 65535 * lut_cam.shape[1] # Corresponding X in the camera image
        img_pts[n][1] = lut_proj[y][x][1] / 65535 * lut_cam.shape[0] # Corresponding Y in the camera image

    assert len(img_pts.shape) == 2 and img_pts.shape[1] == 2


    retval, cam_rotation, cam_translation, inliers = cv2.solvePnPRansac(world_pts, img_pts, cam_internes, cam_disto)

    cam_rotation = cv2.Rodrigues(cam_rotation)[0]

    eprint('# Rotation:')
    eprint_arr(cam_rotation)
    eprint('# Translation:')
    eprint_arr(cam_translation)

    cam_translation = cam_translation.reshape((3,))

    cam_pose = pose_matrix_from_eucl(cam_rotation, cam_translation)

    # Triangulate with the current LUT + merge into lut3d
    pts = triangulate(lut_proj, proj_internes, proj_disto, proj_pose, lut_cam, cam_internes, cam_disto, cam_pose, threshold)
    xy = generate_xy_list(lut_proj, threshold)

    pts[(pts > maxv).any(axis=1)] = (0, 0, 0)
    pts[(pts < minv).any(axis=1)] = (0, 0, 0)
    for n, point in enumerate(pts):
        x, y = xy[n]

        # Replace the current point if the matching distance in the
        # current LUT is smaller than the smallest one to date
        if not (point == 0).all() and confidence[y, x] > lut_proj[y, x, 2]:
            confidence[y, x] = lut_proj[y, x, 2]
            lut3d[y][x] = point

    if debug:
        save_lut3d(i)

# Confidence at each point
if debug:
    confidence[confidence < float('inf')] /= confidence[confidence < float('inf')].max()
    confidence *= 255
    confidence[confidence == float('inf')] = 255

    write_img('confidence.pgm', confidence, 8)


if lut3d_fname:

    h, w = lut3d.shape[0], lut3d.shape[1]

    pts = lut3d.reshape((h * w, 3))

    pts[(pts > maxv).any(axis=1)] = (0, 0, 0)
    pts[(pts < minv).any(axis=1)] = (0, 0, 0)

    out = ((pts - minv)/(maxv - minv)).reshape((h, w, 3)) * 65535

    out[(pts.reshape((h, w, 3)) == 0).all(axis=2)] = (0, 0, 0)

    eprint('min', minv)
    eprint('max', maxv)

    write_img(lut3d_fname, out, 16)

    nonzero_points = lut3d[(lut3d != 0).any(axis=2)]

    assert len(nonzero_points.shape) == 2 and nonzero_points.shape[1] == 3

    np.savetxt(stdout.buffer, nonzero_points)
