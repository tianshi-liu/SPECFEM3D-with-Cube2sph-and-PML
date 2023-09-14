import numpy as np
def xyz2latlond(x, y, z, ASSUME_PERFECT_SPHERE):
    DEGREE_PER_RAD = 180.0 / np.pi
    R = 6371000.0
    r, theta, phi = xyz2rthetaphi(x, y, z)
    theta_out, phi_out = reduce_theta_phi(theta, phi)
    theta_prime = geocentric_2_geographic(theta_out, ASSUME_PERFECT_SPHERE)
    lat = (np.pi / 2 - theta_prime) * DEGREE_PER_RAD
    lon = phi_out * DEGREE_PER_RAD
    depth = R - r
    return lat, lon, depth

def geocentric_2_geographic(theta,ASSUME_PERFECT_SPHERE):
    TINYVAL = 1.0e-9
    FLATTENING_F = 1.0 / 299.8
    ONE_MINUS_F_SQUARED = (1.0 - FLATTENING_F) ** 2
    FACTOR_TAN = 1.0 / ONE_MINUS_F_SQUARED
    if (not ASSUME_PERFECT_SPHERE):
      theta_prime = np.pi / 2 - np.arctan2(FACTOR_TAN * np.cos(theta),  
            np.maximum(TINYVAL, np.sin(theta)))
    else:
      theta_prime = theta
    return theta_prime

def xyz2rthetaphi(x, y, z):
    TINY_VAL = 1e-8
    if (z > - TINY_VAL) and (z <= 0):
        z = - TINY_VAL
    if (z < TINY_VAL) and (z >= 0):
        z = TINY_VAL
    theta = np.arctan2(np.sqrt(x * x + y * y), z)
    if (x > - TINY_VAL) and (x <= 0):
        x = - TINY_VAL
    if (x < TINY_VAL) and (x >= 0):
        x = TINY_VAL
    phi = np.arctan2(y, x)
    r = np.sqrt(x * x + y * y + z * z)
    return r, theta, phi

def reduce_theta_phi(theta, phi):
    TINY_VAL = 1e-8
    if (np.abs(theta) < TINY_VAL):
        theta = theta + TINY_VAL
    if (np.abs(phi) < TINY_VAL):
        phi = phi + TINY_VAL
    if (phi < 0) or (phi > 2 * np.pi):
        cycle = np.abs(int(phi / (2 * np.pi)))
        if (phi < 0):
            phi = phi + (cycle + 1) * 2 * np.pi
        else:
            if (phi > 2 * np.pi):
                phi = phi - cycle * 2 * np.pi
    if (theta < 0) or (theta > np.pi):
        cycle = int(theta / np.pi)
        if (theta > 0):
            if (np.mod(cycle, 2) != 0):
                theta = (cycle + 1) * np.pi - theta
                if (phi < np.pi):
                    phi = phi + np.pi
                else:
                    phi = phi - np.pi
            else:
                theta = theta - cycle * np.pi
        else:
            if (np.mod(cycle, 2) == 0):
                theta = - theta + cycle * np.pi
                if (phi < np.pi):
                    phi = phi + np.pi
                else:
                    phi = phi - np.pi
            else:
                theta = theta - cycle * np.pi
    if (theta < 0) or (theta > np.pi):
        print('theta out of range\n')
    if (phi < 0) or (phi > 2 * np.pi):
        print('phi out of range\n')
    return theta, phi
