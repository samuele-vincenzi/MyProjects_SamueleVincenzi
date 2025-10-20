
import numpy as np


def quaternion_mult(q_rot_1, q_rot_2):
    # The function computes the quaternion multiplication

    # INPUTS:
    # q_rot_1: first quaternion rotation
    # q_rot_2: second quaternion rotation

    # OUTPUT:
    # q: result of the quaternion product

    q1, q2, q3, q4 = q_rot_1
    Q_rot_1 = np.array([
        [q4, -q3, q2, q1],
        [q3, q4, -q1, q2],
        [-q2, q1, q4, q3],
        [-q1, -q2, -q3, q4]
    ])
    q = Q_rot_1 @ q_rot_2
    return q


def quaternion_norm(q):
    # This function normalizes the quaternion

    # INPUT:
    # q: quaternion to be normalized

    # OUTPUT:
    # q_normal: normalized quaternion

    q_normal /= np.linalg.norm(q, axis=0)
    return q_normal
