"""Lidar scans (Hesai OT128, exported as PLY by its viewer): load a capture folder of frames
into one point cloud in the lidar's own frame.

Lidar frame (OT128 manual, 1.3): Z is the rotation axis, up; Y is azimuth 0, the front; X is
to the right of Y (right-handed). Metres, from the lidar's origin (manual 1.5).
A static scan repeats the same scene every frame, so the frames are merged and thinned to
one point per VOXEL_M cube (the brightest)."""
import glob
import os

import numpy as np

VOXEL_M = 0.1
MAX_POINTS = 600_000   # what the browser is sent at most
PLY = np.dtype([("x", "<f4"), ("y", "<f4"), ("z", "<f4"), ("t", "<f8"), ("ring", "<u2"), ("i", "u1")])


def read_ply(path):
    """The vertices of one binary PLY frame (x, y, z, timestamp, ring, intensity)."""
    b = open(path, "rb").read()
    h = b.index(b"end_header\n") + len(b"end_header\n")
    head = b[:h].decode("ascii", "replace").splitlines()
    if "format binary_little_endian 1.0" not in head:
        raise ValueError(f"{os.path.basename(path)}: not a binary little-endian PLY")
    n = int(next(l for l in head if l.startswith("element vertex")).split()[2])
    return np.frombuffer(b, PLY, n, h)


def load_capture(folder, voxel=VOXEL_M):
    """All the PLY frames of a capture folder, merged and thinned: float32 [n, 4] = x, y, z,
    intensity (0-255), in the lidar frame."""
    files = glob.glob(os.path.join(folder, "*.ply"))
    if not files:
        raise ValueError(f"no .ply frames in {folder}")
    a = np.concatenate([read_ply(f) for f in files])
    p = np.stack([a["x"], a["y"], a["z"], a["i"].astype(np.float32)], 1)
    p = p[np.isfinite(p).all(1) & (np.abs(p[:, :3]).sum(1) > 0)]
    p = p[np.argsort(-p[:, 3], kind="stable")]   # brightest first: np.unique keeps the first of each voxel
    _, keep = np.unique(np.floor(p[:, :3] / voxel).astype(np.int64), axis=0, return_index=True)
    p = p[keep]
    if len(p) > MAX_POINTS:
        p = p[np.random.default_rng(0).choice(len(p), MAX_POINTS, replace=False)]
    return p.astype(np.float32), len(files), len(a)

