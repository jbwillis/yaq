import numpy as np
import pytest

from yaq.quaternion import *

def test_skew():
    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    v = np.array([4.0, -2.0, -1.0])

    assert np.allclose(skew(e1) @ e2, e3)
    assert np.allclose(skew(e2) @ e3, e1)
    assert np.allclose(skew(e3) @ e1, e2)
    assert np.allclose(skew(e1).T @ e2, -e3)
    assert np.allclose(skew(e2).T @ e3, -e1)
    assert np.allclose(skew(e3).T @ e1, -e2)
    assert np.allclose(skew(v) @ v, np.zeros(3))

def test_quaternion_product():
    q0 = np.array([1.0, 0.0, 0.0, 0.0])
    assert np.allclose(Lprod(q0), np.eye(4))
    assert np.allclose(Rprod(q0), np.eye(4))

    q1 = np.array([1.0, 1.0, 0.0, 0.0])
    q1 /= np.linalg.norm(q1)
    q2 = np.array([1.0, 1.0, 1.0, 0.0])
    q2 /= np.linalg.norm(q2)
    q3 = np.array([1.0, 1.0, 1.0, 1.0])
    q3 /= np.linalg.norm(q3)

    assert np.allclose(Lprod(q1) @ inv(q1), q0)
    assert np.allclose(Lprod(q2) @ inv(q2), q0)
    assert np.allclose(Lprod(q3) @ inv(q3), q0)

    assert np.allclose(Rprod(q1) @ inv(q1), q0)
    assert np.allclose(Rprod(q2) @ inv(q2), q0)
    assert np.allclose(Rprod(q3) @ inv(q3), q0)

def test_hamilton_convention():
    q1 = np.array([0.0, 1.0, 0.0, 0.0])
    q2 = np.array([0.0, 0.0, 1.0, 0.0])
    q3 = np.array([0.0, 0.0, 0.0, 1.0])

    assert np.allclose(Lprod(q1) @ q2, q3) # ij = k
    assert np.allclose(Lprod(q2) @ q1, -q3) # ji = -k
    assert np.allclose(Lprod(q2) @ q3, q1) # jk = i

    assert np.allclose(Rprod(q1) @ q2, -q3) # ji = -k
    assert np.allclose(Rprod(q2) @ q1, q3) # ij = k
    assert np.allclose(Rprod(q2) @ q3, -q1) # kj = -i


def test_axis_angle():
    q = np.array([1.0, 1.0, 1.0, 1.0])
    q /= np.linalg.norm(q)

    axis, angle = axis_angle_from_q(q)
    qret = q_from_axis_angle(axis, angle)
    assert np.allclose(qret, q)

    for angle in np.deg2rad([-5.0, -195, 45, 90, 135, 270]):
        axis = np.array([1.0, -1.0, 1.0])
        q = q_from_axis_angle(axis, angle)
        ax_ret, an_ret = axis_angle_from_q(q)

        if 0.0 > angle > np.deg2rad(-180.0):
            axis *= -1.0
        elif angle > np.deg2rad(180.0):
            axis *= -1.0

        wrapped = abs(2*np.arctan(np.tan(angle/2)))  # wrap to 0-360deg

        assert np.allclose(ax_ret, axis / np.linalg.norm(axis))
        assert np.allclose(an_ret, wrapped)

def test_axis_angle_identity():
    qret = q_from_axis_angle(np.zeros(3), 0)
    assert np.allclose(qret, np.array([1.0, 0.0, 0.0, 0.0]))
    axis, angle = axis_angle_from_q(np.array([1.0, 0.0, 0.0, 0.0]))
    assert angle == 0

def test_rotm():
    q = np.array([1.0, 0.0, 0.0, 0.0])
    R0 = rotm_from_q(q)
    assert np.allclose(R0, np.eye(3))

    q = np.array([0.0, 1.0, 1.0, 0.0])
    q /= np.linalg.norm(q)

    assert np.allclose(rotm_from_q(q), (Rprod(q).T@Lprod(q))[1:,1:])

    q = np.array([0.0, -1.0, 1.0, 1.0])
    q /= np.linalg.norm(q)

    assert np.allclose(rotm_from_q(q), (Rprod(q).T@Lprod(q))[1:,1:])

def random_unit_q(rng):
    q = rng.normal(size=4)
    return q / np.linalg.norm(q)

_rng = np.random.default_rng(0)
QDOT_CASES = [(random_unit_q(_rng), _rng.normal(size=3)) for _ in range(5)]

@pytest.mark.parametrize("q0, w", QDOT_CASES)
def test_qdot_matches_constant_rate_solution(q0, w):
    # Constant body rate has the exact solution q(t) = q0 ⊗ exp(w t / 2); using a
    # non-identity q0 and random w makes a world-frame (Rprod) qdot fail here.
    def q_at(t):
        return Lprod(q0) @ q_from_axis_angle(w, np.linalg.norm(w) * t)
    h = 1e-6
    q_rate_fd = (q_at(h) - q_at(-h)) / (2 * h)
    assert np.allclose(qdot(q0, w), q_rate_fd, atol=1e-8)

@pytest.mark.parametrize("q0, w", QDOT_CASES)
def test_qdot_gives_body_frame_rotm_rate(q0, w):
    # Ties qdot to the rotation-matrix convention: d/dt(W_R_B) = W_R_B [w_B]x.
    # rotm_from_q is quadratic in q, so the central difference has no truncation error.
    qd = qdot(q0, w)
    h = 1e-4
    rotm_rate_fd = (rotm_from_q(q0 + h * qd) - rotm_from_q(q0 - h * qd)) / (2 * h)
    assert np.allclose(rotm_rate_fd, rotm_from_q(q0) @ skew(w))

@pytest.mark.parametrize("q0, w", QDOT_CASES)
def test_qdot_preserves_norm(q0, w):
    assert np.isclose(q0 @ qdot(q0, w), 0.0)
