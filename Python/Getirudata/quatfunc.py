# Quaternion function definitions

import Ska.Numpy 
from pylab import *
from math import *

def quatmult(q1, q2):
    """Function to multiply two quaternions
       Input q1,q2 : quaternions(5) with q[0] = time
       Output q3 : quaternion product q3 = q1*q2
       q3[0] = 0.0
       q3[1] =  q1[4]*q2[1] - q1[3]*q2[2] + q1[2]*q2[3] + q1[1]*q2[4]
       q3[2] =  q1[3]*q2[1] + q1[4]*q2[2] - q1[1]*q2[3] + q1[2]*q2[4]
       q3[3] = -q1[2]*q2[1] + q1[1]*q2[2] + q1[4]*q2[3] + q1[3]*q2[4] 
       q3[4] = -q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3] + q1[4]*q2[4]
    """
    q3 = zeros(q2.shape)
    q3[0, ] =  q2[0, ]
    q3[1, ] =  q1[4, ] * q2[1, ] - q1[3, ] * q2[2, ] + q1[2, ] * q2[3, ] + q1[1, ] * q2[4, ]
    q3[2, ] =  q1[3, ] * q2[1, ] + q1[4, ] * q2[2, ] - q1[1, ] * q2[3, ] + q1[2, ] * q2[4, ]
    q3[3, ] = -q1[2, ] * q2[1, ] + q1[1, ] * q2[2, ] + q1[4, ] * q2[3, ] + q1[3, ] * q2[4, ] 
    q3[4, ] = -q1[1, ] * q2[1, ] - q1[2, ] * q2[2, ] - q1[3, ] * q2[3, ] + q1[4, ] * q2[4, ]
    return (q3)

def quatconj(q):
    """Function to return quaternion conjugate
       Input q : quaternion(5) with q[0] = time
       Output p : quaternion conjugate
       p[0] =  q[0]
       p[1] = -q[1]
       p[2] = -q[2]
       p[3] = -q[3]
       p[4] =  q[4]
    """
    p = zeros(q.shape)
    p[0, ] =  q[0, ]
    p[1:4, ] = -q[1:4, ]
    p[4, ] =  q[4, ]
    return (p)

def quatnorm(q):
    """Function to return normalized quaternion
       Input q : quaternion(5) with q[0] = time
       Output p : normalized quaternion
    """
    p = zeros(q.shape)
    p[0, ] = q[0, ]
    mag = sqrt(dot(q[1:, ].transpose(), q[1:, ]))
    p[1:, ] = q[1:, ] / mag
    if (p[4, ] < 0.0):
        p[1:, ] = -p[1:, ]
    return (p)

def vect2quat(v):
    """Function to convert a rotation vector to a quaternion
       Input v : vector(4,1) with v[0,] = time
       Output q : quaternion(5) with q[0,] = time
    """
    q = zeros(5, )
    v = v.squeeze()
    q[0, ] = v[0, ]
    vmag = angle = sqrt(dot(v[1:, ].transpose(), v[1:, ]))
    if (angle < 0.0000001): # true if vector conponents small enough for linear approx
        q[1:4, ] = v[1:, ] / 2.0 # only indices 1, 2, 3 being set
        q[4, ] = sqrt(1.0 - dot(q[1:4, ].transpose(), q[1:4, ]))
    else: # true for general case
        q[1:4, ] = v[1:, ] / vmag * sin(angle / 2.0)
        q[4, ] = cos(angle / 2.0)
        if (q[4, ] < 0.0):
            q[1:, ] = -q[1:, ]
    return (q)

def quat2vect(q):
    """Function to convert a quaternion to a rotation vector
       Input q : quaternion(5,1) with q[0,] = time
       Output v : vector(4,1) with v[0,] = time
    """
    v = zeros(q.shape)
    v[0, ] = q[0, ]
    sinang2 = sqrt(dot(q[1:4, ].transpose(), q[1:4, ]))
    angle = 2.0 * atan2(sinang2,q[4, ])
    if (angle < 0.0000001):
        v[1:4, ] = 2.0 * q[1:4, ]
    else:
        v[1:4, ] = q[1:4, ] / sinang2 * angle
    return (v[0:4, ])

def vectmag(v):
    """function to compute magnitude of vector
       Input v : vector of any length
       Output mag : scalar magnitude, ignores 0 component
    """
    mag = sqrt(dot(v[1:, ].transpose(), v[1:, ]))
    return (mag)

def unitvect(v):
    """function to compute unit vector of vector
       Input v : vector of any length
       Output mag : scalar magnitude, ignores 0 component
    """
    u = zeros(shape(v))
    u[0, ] = v[0, ]
    mag = sqrt(dot(v[1:, ].transpose(), v[1:, ]))
    if (mag != 0.0):
        u[1:, ] = v[1:, ] / mag
    return (u)

def quat2mat(q):
    """function to convert quaternion to matrix
       Input quaternion with 0 index = time
       Output matrix(3x3)
    """
    M = zeros((3,3))
    M[0, 0] =  q[1, ] * q[1, ] - q[2, ] * q[2, ] - q[3, ] * q[3, ] + q[4, ] * q[4, ]
    M[0, 1] =  2.0 * (q[1, ] * q[2, ] + q[3, ] * q[4, ])
    M[0, 2] =  2.0 * (q[1, ] * q[3, ] - q[2, ] * q[4, ])
    M[1, 0] =  2.0 * (q[1, ] * q[2, ] - q[3, ] * q[4, ])
    M[1, 1] = -q[1, ] * q[1, ] + q[2, ] * q[2, ] - q[3, ] * q[3, ] + q[4, ] * q[4, ]
    M[1, 2] =  2.0 * (q[2, ] * q[3, ] + q[1, ] * q[4, ])
    M[2, 0] =  2.0 * (q[1, ] * q[3, ] + q[2, ] * q[4, ])
    M[2, 1] =  2.0 * (q[2, ] * q[3, ] - q[1, ] * q[4, ])
    M[2, 2] = -q[1, ] * q[1, ] - q[2, ] * q[2, ] + q[3, ] * q[3, ] + q[4, ] * q[4, ]
    return (M)

def quatxaxis(q):
    """Function to compute the X-axis of an attitude quaternion.
       For a quaternion which transforms from inertial to body coordinates,
       the X-axis is in inertial coordinates.
       Input q : q(5,num), where q(0,num) is time and q(1:,num) is quaternion
       Output X : X(4,num), where X(0,num) is time and X(1:,num) is 3-vector
    """
    X = zeros(q.shape)
    X[0, ] = q[0, ]
    X[1, ] = q[1, ] * q[1, ] - q[2, ] * q[2, ] - q[3, ] * q[3, ] + q[4, ] * q[4, ]
    X[2, ] = 2.0 * (q[1, ] * q[2, ] + q[3, ] * q[4, ])
    X[3, ] = 2.0 * (q[1, ] * q[3, ] - q[2, ] * q[4, ])
    return (X[0:4, ])
