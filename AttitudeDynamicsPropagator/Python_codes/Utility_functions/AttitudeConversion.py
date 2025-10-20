import numpy as np


def attitude_conversion(att, type, q_prev=None):
    # The function computes the conversion among different attitude representations:
    # DCM, Euler angles, quaternion

    # INPUTS:
    # att: attitude vector to be converted
    # type: type of conversion
    # q_prev: previous quaternion to avoid quaternion sign flip

    # OUTPUT:
    # conv_att: converted attitude vector

    match type:

        case 'quat_2_DCM':
            # Conversion from quaternion to DCM
            q1, q2, q3, q4 = att
            conv_att = np.array([
                [q1**2 - q2**2 - q3**2 + q4**2,   2 *
                    (q1*q2 + q3*q4),       2*(q1*q3 - q2*q4)],
                [2*(q1*q2 - q3*q4),              -q1**2 +
                    q2**2 - q3**2 + q4**2, 2*(q2*q3 + q1*q4)],
                [2*(q1*q3 + q2*q4),               2 *
                    (q2*q3 - q1*q4),      -q1**2 - q2**2 + q3**2 + q4**2]
            ])

        case 'DCM_2_quat':
            # Conversion from DCM to quaternion
            DCM = att
            trace = np.trace(DCM)

            if trace > 0:
                q4 = 0.5 * np.sqrt(1 + trace)
                q1 = (DCM[1, 2] - DCM[2, 1]) / (4*q4)
                q2 = (DCM[2, 0] - DCM[0, 2]) / (4*q4)
                q3 = (DCM[0, 1] - DCM[1, 0]) / (4*q4)
            else:
                # Handle cases when trace is small or negative
                if DCM[0, 0] > DCM[1, 1] and DCM[0, 0] > DCM[2, 2]:
                    q1 = 0.5 * np.sqrt(1 + DCM[0, 0] - DCM[1, 1] - DCM[2, 2])
                    q2 = (DCM[0, 1] + DCM[1, 0]) / (4*q1)
                    q3 = (DCM[0, 2] + DCM[2, 0]) / (4*q1)
                    q4 = (DCM[1, 2] - DCM[2, 1]) / (4*q1)
                elif DCM[1, 1] > DCM[2, 2]:
                    q2 = 0.5 * np.sqrt(1 - DCM[0, 0] + DCM[1, 1] - DCM[2, 2])
                    q1 = (DCM[0, 1] + DCM[1, 0]) / (4*q2)
                    q3 = (DCM[1, 2] + DCM[2, 1]) / (4*q2)
                    q4 = (DCM[2, 0] - DCM[0, 2]) / (4*q2)
                else:
                    q3 = 0.5 * np.sqrt(1 - DCM[0, 0] - DCM[1, 1] + DCM[2, 2])
                    q1 = (DCM[0, 2] + DCM[2, 0]) / (4*q3)
                    q2 = (DCM[1, 2] + DCM[2, 1]) / (4*q3)
                    q4 = (DCM[0, 1] - DCM[1, 0]) / (4*q3)

            q = np.array([q1, q2, q3, q4])
            q /= np.linalg.norm(q)  # Normalize

            # Avoid quaternion flip
            if q_prev != None:
                if np.dot(q, q_prev) < 0:
                    q = -q

            conv_att = q

        case 'DCM_2_EulerAngles':
            # Conversion from DCM to Euler angles
            # Rotation 123
            # Singularity at theta = (2n+1)pi/2
            DCM = att
            # Pitch (theta)
            theta = np.arcsin(DCM[2, 0])
            cth = np.cos(theta)
            if abs(cth) > 1e-8:
                psi = np.arctan2(-DCM[1, 0], DCM[0, 0])
                phi = np.arctan2(-DCM[2, 1], DCM[2, 2])
            else:
                # Gimbal lock
                if DCM[2, 0] > 0:   # theta = +pi/2
                    # psi + phi = atan2(C12, -C13)
                    s = np.arctan2(DCM[0, 1], -DCM[0, 2])
                    phi = 0.0
                    psi = s
                else:            # theta = -pi/2
                    # psi - phi = atan2(C12, C13)
                    d = np.arctan2(DCM[0, 1], DCM[0, 2])
                    phi = 0.0
                    psi = d
            conv_att = np.degrees([phi, theta, psi])

        case 'EulerAngles_2_DCM':
            # Conversion from Euler angles to DCM
            phi, theta, psi = att
            phi = np.deg2rad(phi)
            theta = np.deg2rad(theta)
            psi = np.deg2rad(psi)

            c_phi = np.cos(phi)
            s_phi = np.sin(phi)
            c_theta = np.cos(theta)
            s_theta = np.sin(theta)
            c_psi = np.cos(psi)
            s_psi = np.sin(psi)

            DCM = np.array([
                [c_theta * c_psi, c_psi * s_theta * s_phi + s_psi *
                    c_phi, -c_psi * s_theta * c_phi + s_psi * s_phi],
                [-s_psi * c_theta, -s_psi * s_theta * s_phi + c_psi *
                    c_phi, s_psi * s_theta * c_phi + c_psi * s_phi],
                [s_theta, -c_theta * s_phi, c_theta * c_phi],
            ])
            conv_att = DCM

        case 'EulerAngles_2_quat':
            # Conversion from Euler angles to quaternion
            att = attitude_conversion(att, 'EulerAngles_2_DCM')
            conv_att = attitude_conversion(att, 'DCM_2_quat', q_prev)

        case 'quat_2_EulerAngles':
            # Conversion from quaternion to Euler angles
            att = attitude_conversion(att, 'quat_2_DCM')
            conv_att = attitude_conversion(att, 'DCM_2_EulerAngles')

    return conv_att
