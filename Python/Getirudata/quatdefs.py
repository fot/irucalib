# Quaternion function definitions
# from quatdefs import *

import Ska.Numpy 
from pylab import *
from math import *

def vectmag(v):
    """function to compute magnitude of vector
       Input v : vector of any length
       Output mag : scalar magnitude,
    """
    return (np.sqrt((v * v).sum(axis = 0)))

def unitvect(v):
    """Function to compute unit vector of vector array
       Input v : vector of any length
       Output u : array with columns of v unit vectors
    """
    return (v / vectmag(v))

def quatnorm(q):
    """Function to return normalized quaternion
       Input q : quaternion(5) with q[0] = time
       Output q : normalized quaternion
    """
    q = q.copy()
    q = q.reshape(5, -1)
    qmag = vectmag(q[1:, ])
    q[1:, ] = q[1:, ] / qmag
    idx = find(q[4, ] < 0.0)
    q[1:, idx] = -q[1:, idx]
    return (q)

def quatmult(q1, q2):
    """Function to multiply two quaternions
       Input q1,q2 : quaternions(5) with q[0] = time
       Output q3 : quaternion product q3 = q1*q2
       q3[0] = q2[0]
       q3[1] =  q1[4]*q2[1] - q1[3]*q2[2] + q1[2]*q2[3] + q1[1]*q2[4]
       q3[2] =  q1[3]*q2[1] + q1[4]*q2[2] - q1[1]*q2[3] + q1[2]*q2[4]
       q3[3] = -q1[2]*q2[1] + q1[1]*q2[2] + q1[4]*q2[3] + q1[3]*q2[4] 
       q3[4] = -q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3] + q1[4]*q2[4]
    """
    q1 = q1.copy()
    q1 = q1.reshape(5, -1)
    q2 = q2.copy()
    q2 = q2.reshape(5, -1)
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
       q[0] =  q[0]
       q[1] = -q[1]
       q[2] = -q[2]
       q[3] = -q[3]
       q[4] =  q[4]
    """
    q = q.copy()
    q = q.reshape(5, -1)
    q[0, ] =  q[0, ]
    q[1:4, ] = -q[1:4, ]
    q[4, ] =  q[4, ]
    return (q)

def vect2quat(v):
    """Function to convert a rotation vector to a quaternion
       Input v : vector(4,num) with v[0,] = time
       Output q : quaternion(5,num) with q[0,] = time
    """
    v = v.copy()
    v = v.reshape(4, -1)
    q = zeros((5,v.shape[1]))
    q[0, ] = v[0, ]
    vmag = vectmag(v[1:, ])
    idx = find(vmag < 0.0000001)  
    q[1:4, idx] = v[1:4, idx] / 2.0
    q[4, idx] = np.sqrt(1.0 - vmag[idx] * vmag[idx] / 4.0)
    idx = find(vmag >= 0.0000001)  
    q[1:4, idx] = (v[1:4, idx] / vmag[idx]) * np.sin(vmag[idx] / 2.0)
    q[4, idx] = np.cos(vmag[idx] / 2.0)
    return (q)

def quat2vect(q):
    """Function to convert a quaternion to a rotation vector
       Input q : quaternion(5,1) with q[0,] = time
       Output v : vector(4,1) with v[0,] = time
    """
    q = q.copy()
    q = q.reshape(5, -1)
    v = zeros(q.shape)
    v[0, ] = q[0, ]
    sinhalfang = vectmag(q[1:4, ]) # sine of half rotation angle
    angle = 2.0 * np.arctan2(sinhalfang,q[4, ])
    v[1:4, ] = 2.0 * q[1:4, ]
    idx = find(angle >= 0.0000001)
    v[1:4, idx] = q[1:4, idx] / sinhalfang[idx] * angle[idx]
    return (v[0:4, ])

def quat2mat(q):
    """function to convert quaternion to matrix
       Input Q : quaternion(5)
       Output M : matrix(3x3)
    """
    q = q.copy()
    q = q.reshape(5, -1)
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
    q = q.copy()
    q = q.reshape(5, -1)
    X = zeros(q.shape)
    X[0, ] = q[0, ]
    X[1, ] = q[1, ] * q[1, ] - q[2, ] * q[2, ] - q[3, ] * q[3, ] + q[4, ] * q[4, ]
    X[2, ] = 2.0 * (q[1, ] * q[2, ] + q[3, ] * q[4, ])
    X[3, ] = 2.0 * (q[1, ] * q[3, ] - q[2, ] * q[4, ])
    return (X[0:4, ])
    
def crossprod(a,b):
    """Function to compute the cross products of each vector of two vector
       arrays in sequential order.
       Input  a(4,num) : vector array with a[0,:] as time and a[1:,:] as 3-vector
              b(4,num) : vector array with b[0,:] as time and b[1:,:] as 3-vector
       Output c(4,num) : vector array with c[0,:] as time and c[1:,:] as 3-vector
    """
    a = a.copy()
    a = a.reshape(4, -1)
    b = b.copy()
    b = b.reshape(4, -1)
    c = zeros(a.shape)
    c[0, ] = a[0, ]
    c[1, ] = a[2, ] * b[3, ] - a[3, ] * b[2, ]
    c[2, ] = a[3, ] * b[1, ] - a[1, ] * b[3, ]
    c[3, ] = a[1, ] * b[2, ] - a[2, ] * b[1, ]
    return(c)

def quatxaber(q,v):
    """Function to adjust an attitude quaternion for stellar aberration, when
       the stars used to compute the quaternion are not compensated for 
       stellar aberration and are clustered around the X-axis.
       Input  q(5,num): where q[0,:] is time and q[1:,:] is quaternion
              v(4,num): where v[0,:] is time and v[1:,:] is velocity (km/sec)
       Output q(5,num): where q[0,:] is time and q[1:,:] is adjusted quaternion
    """
    q = q.copy()
    q = q.reshape(5, -1)
    v = v.copy()
    v = v.reshape(4, -1)
    sol = 299792.458 # speed of light (km/sec)
    x = quatxaxis(q)
    aberadj = vect2quat(crossprod(x,v/sol))
    q = quatmult(aberadj,q)
    return(q)
    
