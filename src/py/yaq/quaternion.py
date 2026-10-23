import numpy as np

def skew(v):
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]])

def hatq(v):
    return np.block([0.0, v])

def Lprod(q):
    qs = q[0]
    qv = q[1:]
    return np.block([
        [qs, -qv],
        [qv[:,None], qs*np.eye(3) + skew(qv)],
    ])

def Rprod(q):
    qs = q[0]
    qv = q[1:]
    return np.block([
        [qs, -qv],
        [qv[:,None], qs*np.eye(3) - skew(qv)],
    ])

def inv(q):
    qs = q[0]
    qv = q[1:]
    return np.block([qs, -qv])

def rotm_from_q(q):
    qs = q[0]
    qv = q[1:]
    qvx = skew(qv)
    return np.eye(3) + 2*qs*qvx + 2*(qvx@qvx)

def axis_angle_from_q(q):
    if q[0] < 0:
        q = -q

    qs = q[0]
    qv = q[1:]
    qv_norm = np.linalg.norm(qv)

    angle = 2 * np.atan2(qv_norm, qs)

    if qv_norm == 0:
        axis = np.zeros(3)
    else:
        axis = qv / qv_norm

    return axis, angle

def q_from_axis_angle(axis, angle):
    half_angle = angle / 2.0
    axis_norm = np.linalg.norm(axis)
    if not(axis_norm == 0):
        axis = axis / np.linalg.norm(axis)
    qs = np.cos(half_angle)
    qv = np.sin(half_angle) * axis
    return np.block([qs, qv])

def qdot(q, w):
    return 0.5 * Lprod(q) @ hatq(w)