import numpy as np
from scipy.integrate import solve_ivp
from Utility_functions.QuaternionUtilities import quaternion_mult
from Utility_functions.AttitudeConversion import attitude_conversion


def attitude_dynamics_kinematics(t, rot_state, I, M):
    # The function computes the differential equation for the quaternion-angular velocity
    # rotational state

    # INPUTS:
    # t: time instant
    # rot_state: rotational state [quaternion, angular velocity]
    # I: inertia tensor
    # M: external torque

    # OUTPUT:
    # rot_state_dot: time derivative of the rotational state

    q = rot_state[0:4]
    omega = rot_state[4:7]

    # Normalize quaternion to avoid drift
    q /= np.linalg.norm(q)

    # Quaternion kinematics matrix
    omega_4D = np.hstack([omega, 0.0])
    quat_product = quaternion_mult(q, omega_4D)

    # Quaternion derivative
    q_dot = 0.5 * quat_product

    # Rotational dynamics (Euler's equation)
    I_inv = np.linalg.inv(I)
    omega_dot = I_inv @ (M - np.cross(omega, I @ omega))

    # Concatenate derivatives
    rot_state_dot = np.concatenate((q_dot, omega_dot))
    return rot_state_dot


def attitude_propagator(I, M, tvec, state_in, method, att_type_in):
    # The function propagates the rotational state

    # INPUTS:
    # I: inertia tensor
    # M: external torque
    # tvec: vector containing the time grid
    # state_in: initial rotational state [quaternion, angular velocity]
    # method: numerical integration method for 'solve_ivp'
    # att_type_in: initial attitude representation (DCM, Euler angles or quaternion)

    # OUTPUT:
    # q_vec: time evolution of the quaternion
    # eul_angles: time evolution of the euler angles
    # DCM: time evolution of the DCM
    # omega_vec: time evolution of the angular rates

    # Convert initial attitude representation to quaternion
    match att_type_in:

        case 'DCM':
            att_DCM = state_in[0:10]
            att_quat = attitude_conversion(att_DCM, 'DCM_2_quat', q_prev=None)
            state_in = np.concatenate([att_quat, state_in[10:12]])

        case 'Euler angles':
            att_eul = state_in[0:10]
            att_quat = attitude_conversion(
                att_eul, 'EulerAngles_2_quat', q_prev=None)
            state_in = np.concatenate([att_quat, state_in[3:6]])

        case 'quaternion':
            state_in = state_in

    # Propagate the rotational state
    t_span = (tvec[0], tvec[-1])
    args = (I, M)
    info = solve_ivp(
        attitude_dynamics_kinematics,
        t_span,
        state_in,
        args=args,
        t_eval=tvec,
        method=method,
        rtol=1e-9,
        atol=1e-12,
    )
    q_vec = info.y[0:4, :]
    q_vec /= np.linalg.norm(q_vec, axis=0)
    omega_vec = info.y[4:7, :]

    # Compute euler angles and DCM
    eul_angles = np.zeros([3, len(q_vec[0, :])])
    DCM = np.zeros([3, 3, len(q_vec[0, :])])
    for ii in range(len(q_vec[0, :])):
        if ii == 0:
            eul_angles[:, ii] = attitude_conversion(
                q_vec[:, ii], 'quat_2_EulerAngles')
            DCM[:, :, ii] = attitude_conversion(
                q_vec[:, ii], 'quat_2_DCM')
        else:
            eul_angles[:, ii] = attitude_conversion(
                q_vec[:, ii], 'quat_2_EulerAngles', q_prev=q_vec[:, ii-1])
            DCM[:, :, ii] = attitude_conversion(
                q_vec[:, ii], 'quat_2_DCM', q_prev=q_vec[:, ii-1])
    eul_angles = np.unwrap(eul_angles, axis=1)

    return q_vec, eul_angles, DCM, omega_vec
