from __future__ import print_function
import numpy as np
import itertools as it

def inconsistent((g1, g2)):
    return (g1-g2)**2 == 4

def find_agreement(x1, x2):
    w = 1000
    seg_1 = x1[:w]
    seg_2 = x2[:w]
    disagreements = np.where(map(inconsistent, it.izip(seg_1, seg_2)))[0]
    while len(disagreements) == 0 and w < len(x1):
        w*=2
        seg_1 = x1[:w]
        seg_2 = x2[:w]
        disagreements = np.where(map(inconsistent, it.izip(seg_1, seg_2)))[0]
    if len(disagreements) > 0:
        return disagreements[0]
    else:
        return len(x1)

def find_boundaries(x1, x2, i):
    right_boundary = find_agreement(x1[i:], x2[i:]) + i
    left_boundary = i - find_agreement(x1[i::-1], x2[i::-1]) + 1
