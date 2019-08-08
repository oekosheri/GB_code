#!/usr/bin/env python

"""
This module is a collection of functions that produce CSL properties.
When run from the terminal, the code runs in two modes.

 First mode:
  'csl_generator.py u v w [limit]' ----->  Where the u v w are the
  indices of the rotation axis such as 1 0 0, 1 1 1, 1 1 0 and so on. The limit
  is the maximum Sigma of interest.
  (the limit by default: 100)

 Second mode:
 'csl_generator.py u v w basis sigma [limit]' -----> Where basis is
  either fcc, bcc, diamond or sc. You read the sigma of interest from the first
  mode run. The limit here refers to CSL GB inclinations. The bigger the limit,
  the higher the indices of CSL planes.
  (the limit by default: 2)

  your chosen axis, basis and sigma will be written to an io_file which will
  later be read by gb_geberator.py.
"""

import sys
import random
from math import degrees, atan, sqrt, pi, ceil, cos, acos, sin, gcd, radians
import numpy as np
from numpy import dot, cross
from numpy.linalg import det, norm, inv


def get_cubic_sigma(uvw, m, n=1):
    """
    CSL analytical formula based on the book:
    'Interfaces in crystalline materials',
     Sutton and Balluffi, clarendon press, 1996.
     generates possible sigma values.
    arguments:
    uvw -- the axis
    m,n -- two integers (n by default 1)

    """
    u, v, w = uvw
    sqsum = u*u + v*v + w*w
    sigma = m*m + n*n * sqsum
    while sigma != 0 and sigma % 2 == 0:
        sigma /= 2
    return sigma if sigma > 1 else None


def get_cubic_theta(uvw, m, n=1):
    """
    generates possible theta values.
    arguments:
    uvw -- the axis
    m,n -- two integers (n by default 1)
    """
    u, v, w = uvw
    sqsum = u*u + v*v + w*w
    if m > 0:
        return 2 * atan(sqrt(sqsum) * n / m)
    else:
        return pi


def get_theta_m_n_list(uvw, sigma):
    """
    Finds integers m and n lists that match the input sigma.
    """
    if sigma == 1:
        return [(0., 0., 0.)]
    thetas = []
    max_m = int(ceil(sqrt(4*sigma)))

    for m in range(1, max_m):
        for n in range(1, max_m):
            if gcd(m, n) == 1:
                s = get_cubic_sigma(uvw, m, n)
            if s == sigma:
                theta = (get_cubic_theta(uvw, m, n))
                thetas.append((theta, m, n))
                thetas.sort(key=lambda x: x[0])
    return thetas


def print_list(uvw, limit):
    """
    prints a list of smallest sigmas/angles for a given axis(uvw).
    """
    for i in range(limit):
        tt = get_theta_m_n_list(uvw, i)
        if len(tt) > 0:
            theta, _, _ = tt[0]
            print("Sigma:   {0:3d}  Theta:  {1:5.2f} "
                  .format(i, degrees(theta)))


def rot(a, Theta):
    """
    produces a rotation matrix.
    arguments:
    a -- an axis
    Theta -- an angle
    """
    c = cos(Theta)
    s = sin(Theta)
    a = a / norm(a)
    ax, ay, az = a
    return np.array([[c + ax * ax * (1 - c), ax * ay * (1 - c) - az * s,
                      ax * az * (1 - c) + ay * s],
                    [ay * ax * (1 - c) + az * s, c + ay * ay * (1 - c),
                        ay * az * (1 - c) - ax * s],
                     [az * ax * (1 - c) - ay * s, az * ay * (1 - c) + ax * s,
                      c + az * az * (1 - c)]])


# Helpful Functions:
# -------------------#

def integer_array(A, tol=1e-7):
    """
    returns True if an array is ineteger.
    """
    return np.all(abs(np.round(A) - A) < tol)


def angv(a, b):
    """
    returns the angle between two vectors.
    """
    ang = acos(np.round(dot(a, b)/norm(a)/norm(b), 8))
    return round(degrees(ang), 7)


def ang(a, b):
    """
    returns the cos(angle) between two vectors.
    """
    ang = np.round(dot(a, b)/norm(a)/norm(b), 7)
    return abs(ang)


def CommonDivisor(a):
    """
    returns the common factor of vector a and the reduced vector.
    """
    CommFac = []
    a = np.array(a)
    for i in range(2, 100):
        while (a[0] % i == 0 and a[1] % i == 0 and a[2] % i == 0):
            a = a / i
            CommFac.append(i)
    return(a.astype(int), (abs(np.prod(CommFac))))


def SmallestInteger(a):
    """
    returns the smallest multiple integer to make an integer array.
    """
    a = np.array(a)
    for i in range(1, 200):
        testV = i * a
        if integer_array(testV):
            break
    return (testV, i) if integer_array(testV) else None


def integerMatrix(a):
    """
    returns an integer matrix from row vectors.
    """
    Found = True
    b = np.zeros((3, 3))
    a = np.array(a)
    for i in range(3):
        for j in range(1, 2000):
            testV = j * a[i]
            if integer_array(testV):
                b[i] = testV
                break
        if all(b[i] == 0):
            Found = False
            print("Can not make integer matrix!")
    return (b) if Found else None


def SymmEquivalent(arr):
    """
    returns cubic symmetric eqivalents of the given 2 dimensional vector.
    """
    Sym = np.zeros([24, 3, 3])
    Sym[0, :] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Sym[1, :] = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    Sym[2, :] = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
    Sym[3, :] = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
    Sym[4, :] = [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]
    Sym[5, :] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    Sym[6, :] = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
    Sym[7, :] = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
    Sym[8, :] = [[-1, 0, 0], [0, 0, -1], [0, -1, 0]]
    Sym[9, :] = [[-1, 0, 0], [0, 0, 1], [0, 1, 0]]
    Sym[10, :] = [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
    Sym[11, :] = [[1, 0, 0], [0, 0, 1], [0, -1, 0]]
    Sym[12, :] = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    Sym[13, :] = [[0, 1, 0], [0, 0, -1], [-1, 0, 0]]
    Sym[14, :] = [[0, -1, 0], [0, 0, 1], [-1, 0, 0]]
    Sym[15, :] = [[0, -1, 0], [0, 0, -1], [1, 0, 0]]
    Sym[16, :] = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    Sym[17, :] = [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]
    Sym[18, :] = [[0, 0, -1], [1, 0, 0], [0, -1, 0]]
    Sym[19, :] = [[0, 0, -1], [-1, 0, 0], [0, 1, 0]]
    Sym[20, :] = [[0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    Sym[21, :] = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
    Sym[22, :] = [[0, 0, 1], [0, -1, 0], [1, 0, 0]]
    Sym[23, :] = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]

    arr = np.atleast_2d(arr)
    Result = []
    for i in range(len(Sym)):
        for j in range(len(arr)):
            Result.append(dot(Sym[i, :], arr[j]))
    Result = np.array(Result)
    return np.unique(Result, axis=0)


def Tilt_Twist_comp(v1, uvw, m, n):
    """
    returns the tilt and twist components of a given GB plane.
    arguments:
    v1 -- given gb plane
    uvw -- axis of rotation
    m,n -- the two necessary integers
    """
    theta = get_cubic_theta(uvw, m, n)
    R = rot(uvw, theta)
    v2 = np.round(dot(R, v1), 6).astype(int)
    tilt = angv(v1, v2)
    if abs(tilt - degrees(theta)) < 10e-5:
        print("Pure tilt boundary with a tilt component: {0:6.2f}"
              .format(tilt))
    else:
        twist = 2 * acos(cos(theta / 2) / cos(radians(tilt / 2)))
        print("Tilt component: {0:<6.2f} Twist component: {1:6.2f}"
              .format(tilt, twist))


def Create_Possible_GB_Plane_List(uvw, m=5, n=1, lim=5):
    """
    generates GB planes and specifies the character.

    arguments:
    uvw -- axis of rotation.
    m,n -- the two necessary integers
    lim -- upper limit for the plane indices

    """
    uvw = np.array(uvw)
    Theta = get_cubic_theta(uvw, m, n)
    Sigma = get_cubic_sigma(uvw, m, n)
    R1 = rot(uvw, Theta)

    # List and character of possible GB planes:
    x = np.arange(-lim, lim + 1, 1)
    y = x
    z = x
    indice = (np.stack(np.meshgrid(x, y, z)).T).reshape(len(x) ** 3, 3)
    indice_0 = indice[np.where(np.sum(abs(indice), axis=1) != 0)]
    indice_0 = indice_0[np.argsort(norm(indice_0, axis=1))]

    # extract the minimal cell:
    Min_1, Min_2 = Create_minimal_cell_Method_1(Sigma, uvw, R1)
    V1 = np.zeros([len(indice_0), 3])
    V2 = np.zeros([len(indice_0), 3])
    GBtype = []
    tol = 0.001
    # Mirrorplanes cubic symmetry
    MP = np.array([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1],
                   [1, 1, 0],
                   [0, 1, 1],
                   [1, 0, 1],
                   ], dtype='float')
    # Find GB plane coordinates:
    for i in range(len(indice_0)):
        if CommonDivisor(indice_0[i])[1] <= 1:
            V1[i, :] = (indice_0[i, 0] * Min_1[:, 0] +
                        indice_0[i, 1] * Min_1[:, 1] +
                        indice_0[i, 2] * Min_1[:, 2])
            V2[i, :] = (indice_0[i, 0] * Min_2[:, 0] +
                        indice_0[i, 1] * Min_2[:, 1] +
                        indice_0[i, 2] * Min_2[:, 2])

    V1 = (V1[~np.all(V1 == 0, axis=1)]).astype(int)
    V2 = (V2[~np.all(V2 == 0, axis=1)]).astype(int)
    MeanPlanes = (V1 + V2) / 2

    # Check the type of GB plane: Symmetric tilt, tilt, twist
    for i in range(len(V1)):
        if ang(V1[i], uvw) < tol:

            for j in range(len(SymmEquivalent(MP))):
                if 1 - ang(MeanPlanes[i], SymmEquivalent(MP)[j]) < tol:
                    GBtype.append('Symmetric Tilt')
                    break
            else:
                GBtype.append('Tilt')
        elif 1 - ang(V1[i], uvw) < tol:
            GBtype.append('Twist')
        else:
            GBtype.append('Mixed')

    return (V1, V2, MeanPlanes, GBtype)


def Create_minimal_cell_Method_1(sigma, uvw, R):
    """
    finds Minimal cell by means of a numerical search.
    (An alternative analytical method can be used too).
    arguments:
    sigma -- gb sigma
    uvw -- rotation axis
    R -- rotation matrix
    """
    uvw = np.array(uvw)
    MiniCell_1 = np.zeros([3, 3])
    MiniCell_1[:, 2] = uvw

    lim = 20
    x = np.arange(-lim, lim + 1, 1)
    y = x
    z = x
    indice = (np.stack(np.meshgrid(x, y, z)).T).reshape(len(x) ** 3, 3)

    # remove 0 vectors and uvw from the list
    indice_0 = indice[np.where(np.sum(abs(indice), axis=1) != 0)]
    condition1 = ((abs(dot(indice_0, uvw) / norm(indice_0, axis=1) /
                       norm(uvw))).round(7))
    indice_0 = indice_0[np.where(condition1 != 1)]

    if MiniCell_search(indice_0, MiniCell_1, R, sigma):

        M1, M2 = MiniCell_search(indice_0, MiniCell_1, R, sigma)
        return (M1, M2)
    else:
        return None


def MiniCell_search(indices, MiniCell_1, R, sigma):

    tol = 0.001
    # norm1 = norm(indices, axis=1)
    newindices = dot(R, indices.T).T
    nn = indices[np.all(abs(np.round(newindices) - newindices) < 1e-6, axis=1)]
    TestVecs = nn[np.argsort(norm(nn, axis=1))]
    # print(len(indices), len(TestVecs),TestVecs[:20])

    Found = False
    count = 0
    while (not Found) and count < len(TestVecs) - 1:
        if 1 - ang(TestVecs[count], MiniCell_1[:, 2]) > tol:
            # and  (ang(TestVecs[i],uvw) > tol):
            MiniCell_1[:, 1] = (TestVecs[count])
            count += 1
            for j in range(len(TestVecs)):
                if (1 - ang(TestVecs[j], MiniCell_1[:, 2]) > tol) and (
                        1 - ang(TestVecs[j], MiniCell_1[:, 1]) > tol):
                    if (ang(TestVecs[j],
                            cross(MiniCell_1[:, 2], MiniCell_1[:, 1])) > tol):
                        #  The condition that the third vector can not be
                        #  normal to any other two.
                        #  and (ang(TestVecs[i],uvw)> tol) and
                        # (ang(TestVecs[i],MiniCell[:,1])> tol)):
                        MiniCell_1[:, 0] = (TestVecs[j]).astype(int)
                        Det1 = abs(round(det(MiniCell_1), 5))
                        MiniCell_1 = (MiniCell_1).astype(int)
                        MiniCell_2 = ((np.round(dot(R, MiniCell_1), 7))
                                      .astype(int))
                        Det2 = abs(round(det(MiniCell_2), 5))

                        if ((abs(Det1 - sigma)) < tol and
                                (abs(Det2 - sigma)) < tol):
                            Found = True
                            break
    if Found:
        return MiniCell_1, MiniCell_2
    else:
        return Found


def Basis(basis):
    """
    defines the basis.
    """
    # Cubic basis
    if str(basis) == 'fcc':
        basis = np.array([[0, 0, 0],
                          [0.5, 0, 0.5],
                          [0.5, 0.5, 0],
                          [0, 0.5, 0.5]], dtype=float)
    elif str(basis) == 'bcc':
        basis = np.array([[0, 0, 0],
                          [0.5, 0.5, 0.5]], dtype=float)
    elif str(basis) == 'sc':
        basis = np.eye(3)

    elif str(basis) == 'diamond':
        basis = np.array([[0, 0, 0],
                          [0.5, 0, 0.5],
                          [0.5, 0.5, 0],
                          [0, 0.5, 0.5],
                          [0.25, 0.25, 0.25],
                          [0.75, 0.25, 0.75],
                          [0.75, 0.75, 0.25],
                          [0.25, 0.75, 0.75]], dtype=float)
    else:
        print('Sorry! For now only works for cubic lattices ...')
        sys.exit()

    return basis


def Find_Orthogonal_cell(basis, uvw, m, n, GB1):
    """
    finds Orthogonal cells from the CSL minimal cells.
    arguments:
    basis -- lattice basis
    uvw -- rotation axis
    m,n -- two necessary integers
    GB1 -- input plane orientation
    """
    # Changeable limit
    lim = 15
    uvw = np.array(uvw)
    Theta = get_cubic_theta(uvw, m, n)
    Sigma = get_cubic_sigma(uvw, m, n)
    R = rot(uvw, Theta)
    GB2 = np.round((dot(R, GB1)), 6)
    x = np.arange(-lim, lim + 1, 1)
    y = x
    z = x
    indice = (np.stack(np.meshgrid(x, y, z)).T).reshape(len(x) ** 3, 3)
    indice_0 = indice[np.where(np.sum(abs(indice), axis=1) != 0)]
    indice_0 = indice_0[np.argsort(norm(indice_0, axis=1))]
    OrthoCell_1 = np.zeros([3, 3])
    OrthoCell_1[:, 0] = np.array(GB1)
    OrthoCell_2 = np.zeros([3, 3])
    OrthoCell_2[:, 0] = np.array(GB2)
    # extract the minimal cells:
    Min_1, Min_2 = Create_minimal_cell_Method_1(Sigma, uvw, R)
    # Find Ortho vectors:
    tol = 0.001
    if ang(OrthoCell_1[:, 0], uvw) < tol:
        OrthoCell_1[:, 1] = uvw
        OrthoCell_2[:, 1] = uvw
    else:

        for i in range(len(indice_0)):

            v1 = (indice_0[i, 0] * Min_1[:, 0] +
                  indice_0[i, 1] * Min_1[:, 1] +
                  indice_0[i, 2] * Min_1[:, 2])
            v2 = (indice_0[i, 0] * Min_2[:, 0] +
                  indice_0[i, 1] * Min_2[:, 1] +
                  indice_0[i, 2] * Min_2[:, 2])
            if ang(v1, OrthoCell_1[:, 0]) < tol:
                OrthoCell_1[:, 1] = v1
                OrthoCell_2[:, 1] = v2
                break
    OrthoCell_1[:, 2] = np.cross(OrthoCell_1[:, 0], OrthoCell_1[:, 1])
    OrthoCell_2[:, 2] = np.cross(OrthoCell_2[:, 0], OrthoCell_2[:, 1])

    if (CommonDivisor(OrthoCell_1[:, 2])[1] ==
            CommonDivisor(OrthoCell_2[:, 2])[1]):
        OrthoCell_1[:, 2] = CommonDivisor(OrthoCell_1[:, 2])[0]
        OrthoCell_2[:, 2] = CommonDivisor(OrthoCell_2[:, 2])[0]
    # OrthoCell_1 = OrthoCell_1
    # OrthoCell_2 = OrthoCell_2
    # Test
    # OrthoCell_3 = (dot(R, OrthoCell_1))
    Volume_1 = (round(det(OrthoCell_1), 5))
    Volume_2 = (round(det(OrthoCell_2), 5))
    Num = Volume_1 * len(Basis(basis)) * 2

    if Volume_1 == Volume_2:
        OrthoCell_1 = OrthoCell_1.astype(float)
        OrthoCell_2 = OrthoCell_2.astype(float)

        if basis == 'sc' or basis == 'diamond':

            return ((OrthoCell_1.astype(float),
                     OrthoCell_2.astype(float), Num.astype(int)))

        elif basis == 'fcc' or basis == 'bcc':
            ortho1, ortho2 = Ortho_fcc_bcc(basis, OrthoCell_1, OrthoCell_2)
            Volume_1 = (round(det(ortho1), 5))
            Num = Volume_1 * len(Basis(basis)) * 2
            return (ortho1, ortho2, Num.astype(int))

    else:
        return None


def print_list_GB_Planes(uvw, basis, m, n, lim=3):
    """
    prints lists of GB planes given an axis, basis, m and n.
    """
    uvw = np.array(uvw)
    V1, V2, _, Type = Create_Possible_GB_Plane_List(uvw, m, n, lim)
    for i in range(len(V1)):
        Or = Find_Orthogonal_cell(basis, uvw, m, n, V1[i])
        if Or:
            print("{0:<20s}   {1:<20s}   {2:<20s}   {3:<10s}"
                  .format(str(V1[i]), str(V2[i]), Type[i], str(Or[2])))


# ___CSL/DSC vector construction___#

# According to Grimmer et al. (Acta Cryst. (1974). A30, 197-207) ,
#  DSC and CSL lattices for bcc and fcc were made from the Sc lattice
# via body_centering and face_centering:
def odd_even(M1):
    """
    finds odd and even elements of a matrix.
    """
    d_e = np.array([['a', 'a', 'a'],
                    ['a', 'a', 'a'],
                    ['a', 'a', 'a']], dtype=str)
    for i in range(3):
        for j in range(3):
            if abs(M1[i][j]) % 2 == 0:
                d_e[i][j] = 'e'
            else:
                d_e[i][j] = 'd'
    return d_e


def self_test_b(a):
    z_b = np.eye(3, 3)
    M = a.copy()
    for i in range(3):
        if np.all(odd_even(M)[:, i] == ['d', 'd', 'd']):
            z_b[i][i] = 0.5
            break

    return z_b


def binary_test_b(a):
    count = 0
    z_b = np.eye(3, 3)
    for i in [0, 1]:
        for j in [1, 2]:
            if i != j and count < 1:
                M = a.copy()
                M[:, j] = M[:, i] + M[:, j]
                if np.all(odd_even(M)[:, j] == ['d', 'd', 'd']):
                    count = count + 1
                    z_b[i][j] = 0.5
                    z_b[j][j] = 0.5
    return z_b


def tertiary_test_b(a):
    z_b = np.eye(3, 3)
    M = a.copy()
    M[:, 2] = M[:, 0] + M[:, 1] + M[:, 2]
    for k in range(3):
        if np.all(odd_even(M)[:, k] == ['d', 'd', 'd']):
            z_b[0][k] = 0.5
            z_b[1][k] = 0.5
            z_b[2][k] = 0.5
            break
    return z_b


def body_centering(b):
    """
    converting a single crystal minimal cell to a bcc one.
    """
    z_b = np.eye(3, 3)
    while det(z_b) != 0.5:
        if det(self_test_b(b)) == 0.5:
            z_b = self_test_b(b)
            break
        if det(binary_test_b(b)) == 0.5:
            z_b = binary_test_b(b)
            break
        if det(tertiary_test_b(b)) == 0.5:
            z_b = tertiary_test_b(b)
            break
    return z_b


def face_centering(a):
    """
    converting a single crystal minimal cell to an fcc one.
    """
    z_f = np.eye(3, 3)
    M = a.copy()
    count = 0
    for i in range(3):
        if (np.all(odd_even(M)[:, i] == ['d', 'd', 'e']) or
                np.all(odd_even(M)[:, i] == ['e', 'd', 'd']) or
                np.all(odd_even(M)[:, i] == ['d', 'e', 'd'])):
            count = count + 1
            z_f[i][i] = 0.5
    if det(z_f) == 0.25:
        return (z_f)
    else:
        for i in [0, 1]:
            for j in [1, 2]:
                if i != j and count < 2:
                    M = a.copy()
                    M[:, j] = M[:, i] + M[:, j]
                    if (np.all(odd_even(M)[:, j] == ['d', 'd', 'e']) or
                            np.all(odd_even(M)[:, j] == ['e', 'd', 'd']) or
                            np.all(odd_even(M)[:, j] == ['d', 'e', 'd'])):
                        count = count + 1
                        z_f[i][j] = 0.5
                        z_f[j][j] = 0.5

    return z_f if det(z_f) == 0.25 else None


def DSC_vec(basis, sigma, minicell):
    """
    a discrete shift complete (DSC)
    network for the given sigma and minimal cell.
    arguments:
    basis -- a lattice basis(fcc or bcc)
    sigma -- gb sigma
    minicell -- gb minimal cell
    """
    D_sc = np.round(sigma * (inv(minicell).T), 6).astype(int)
    if basis == 'sc':
        D = D_sc.copy()
    if basis == 'bcc':
        D = dot(D_sc, body_centering(D_sc))
    if basis == 'fcc' or basis == 'diamond':
        D = dot(D_sc, face_centering(D_sc))
    return D


def CSL_vec(basis, minicell):
    """
    CSL minimal cell for sc, fcc and bcc.
    arguments:
    basis -- a lattice basis(fcc or bcc)
    minicell -- gb minimal cell
    """
    C_sc = minicell.copy()
    if basis == 'sc':
        C = C_sc.copy()
    if basis == 'bcc':
        C = dot(C_sc, body_centering(C_sc))
    if basis == 'fcc':
        C = dot(C_sc, face_centering(C_sc))
    return C


def DSC_on_plane(D, Pnormal):
    """
    projects the given DSC network on a given plane.
    """
    D_proj = np.zeros((3, 3))
    Pnormal = np.array(Pnormal)
    for i in range(3):
        D_proj[:, i] = (D[:, i] - (dot(D[:, i], Pnormal) / norm(Pnormal)) *
                        Pnormal/norm(Pnormal))
    return D_proj


def CSL_density(basis, minicell, plane):
    """
    returns the CSL density of a given plane and its d_spacing.
    """
    plane = np.array(plane)
    C = CSL_vec(basis, minicell)
    h = dot(C.T, plane)
    h = SmallestInteger(h)[0]
    h = CommonDivisor(h)[0]
    G = inv(dot(C.T, C))
    hnorm = np.sqrt(dot(h.T, dot(G, h)))
    density = 1/(hnorm * det(C))
    return (abs(density), 1/hnorm)

# An auxilary function to help reduce the size of the small orthogonal cell
# that is decided otherwise based on Sc for fcc and bcc


def Ortho_fcc_bcc(basis, O1, O2):
    ortho1 = np.zeros((3, 3))
    ortho2 = np.zeros((3, 3))
    if basis == 'fcc':
        base = np.delete(Basis('fcc').T, 0, 1)
    elif basis == 'bcc':
        base = ((np.array([[0.5, 0.5, 0.5], [0.5, 0.5, -0.5],
                           [-0.5, 0.5, 0.5]])).T)

    for i in range(3):
        Min_d = min(CommonDivisor(dot(O1[:, i], inv(base)))[1],
                    CommonDivisor(dot(O2[:, i], inv(base)))[1])
        ortho1[:, i] = O1[:, i] / Min_d
        ortho2[:, i] = O2[:, i] / Min_d
    return (ortho1, ortho2)
# Writing to a yaml file that will be read by gb_generator


def Write_to_io(axis, m, n, basis):
    """
    an input file for gb_generator.py that can be customized.
    It also contains the output from the usage of csl_generator.py.
    """

    my_dict = {'GB_plane': str([axis[0], axis[1], axis[2]]),
               'lattice_parameter': '4',
               'overlap_distance': '0.0', 'which_g': 'g1',
               'rigid_trans': 'no', 'a': '10', 'b': '5',
               'dimensions': '[1,1,1]',
               'File_type': 'LAMMPS'}

    with open('io_file', 'w') as f:
        f.write('### input parameters for gb_generator.py ### \n')
        f.write('# CSL plane of interest that you read from the output of '
                'csl_generator as GB1 \n')
        f.write(list(my_dict.keys())[0] + ': ' + list(my_dict.values())[0] +
                '\n\n')
        f.write('# lattice parameter in Angstrom \n')
        f.write(list(my_dict.keys())[1] + ': ' + list(my_dict.values())[1] +
                '\n\n')
        f.write('# atoms that are closer than this fraction of the lattice '
                'parameter will be removed \n')
        f.write('# either from grain1 (g1) or from grain2 (g2). If you choose '
                '0 no atoms will be removed \n')
        f.write(list(my_dict.keys())[2] + ': ' + list(my_dict.values())[2] +
                '\n\n')
        f.write('# decide which grain the atoms should be removed from \n')
        f.write(list(my_dict.keys())[3]+': ' + str(list(my_dict.values())[3]) +
                '\n\n')
        f.write('# decide whether you want rigid body translations to be done '
                'on the GB_plane or not (yes or no)\n')

        f.write('# When yes, for any GB aside from twist GBs, the two inplane \n'
        '# CSL vectors will be divided by integers a and b to produce a*b initial \n'
        '# configurations. The default values produce 50 initial structures \n'
        '# if you choose no for rigid_trans, you do not need to care about a and b. \n'
        '# twist boundaries are handled internally \n')

        f.write(list(my_dict.keys())[4] + ': ' +
                str(list(my_dict.values())[4]) + '\n')

        f.write(list(my_dict.keys())[5] + ': ' +
                str(list(my_dict.values())[5]) + '\n')

        f.write(list(my_dict.keys())[6] + ': ' +
                str(list(my_dict.values())[6]) + '\n\n')

        f.write('# dimensions of the supercell in: [l1,l2,l3],  where l1 is'
                'the direction along the GB_plane normal\n')
        f.write('#  and l2 and l3 are inplane dimensions \n')
        f.write(list(my_dict.keys())[7] + ': ' + list(my_dict.values())[7] +
                '\n\n')
        f.write('# File type, either VASP or LAMMPS input \n')
        f.write(list(my_dict.keys())[8] + ': ' + list(my_dict.values())[8] +
                '\n\n\n')
        f.write('# The following is your csl_generator output.'
                ' YOU DO NOT NEED TO CHANGE THEM! \n\n')
        f.write('axis'+': ' + str([axis[0], axis[1], axis[2]]) + '\n')
        f.write('m' + ': ' + str(m) + '\n')
        f.write('n' + ': ' + str(n) + '\n')
        f.write('basis' + ': ' + str(basis) + '\n')

    f.close()
    return


def main():

    if (len(sys.argv) != 4 and len(sys.argv) != 5 and len(sys.argv) != 6 and
            len(sys.argv) != 7):
        print(__doc__)
    else:
        uvw = np.array([int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])])
        uvw = CommonDivisor(uvw)[0]

    if len(sys.argv) == 4:
        limit = 100
        print("   List of possible CSLs for {} axis sorted by Sigma   "
              .format(str(uvw)))
        print_list(uvw, limit)
        print("\n Choose a basis, pick a sigma and use the second mode!\n")

    if len(sys.argv) == 5:

        try:
            limit = int(sys.argv[4])
            print("    List of possible CSLs for {} axis sorted by Sigma   "
                  .format(str(uvw)))
            print_list(uvw, limit)
            print("\n Choose a basis, pick a sigma and use the second mode!\n")

        except:
            print("""
                Careful! limit must be an integer ...
                  """, __doc__)

    if len(sys.argv) == 6:

        lim = 2
        basis = str(sys.argv[4])
        sigma = int(sys.argv[5])

        try:
            _, m, n = get_theta_m_n_list(uvw, sigma)[0]
            Write_to_io(uvw, m, n, basis)

            print("----------List of possible CSL planes for Sigma {}---------"
                  .format(sigma))
            print(" GB1-------------------GB2-------------------Type----------"
                  "Number of Atoms ")
            print_list_GB_Planes(uvw, basis, m, n, lim)
            print(" \nPick a GB plane and customize the io_file! ")
            print(" then run : python gb_generator.py io_file\n ")

        except:
            print("Your input sigma is wrong!")
            sys.exit()

    if len(sys.argv) == 7:

        basis = str(sys.argv[4])
        sigma = int(sys.argv[5])

        try:

            _, m, n = get_theta_m_n_list(uvw, sigma)[0]
            lim = int(sys.argv[6])

            if lim > 10:
                print(2*'\n')
                print('You have chosen a large limit! It may take a while ...')
                print(2*'\n')

            print("----------List of possible CSL planes for Sigma {}---------"
                  .format(sigma))
            print(" GB1-------------------GB2-------------------Type----------"
                  "Number of Atoms ")

            print_list_GB_Planes(uvw, basis, m, n, lim)
            print(" \nPick a GB plane and customize the io_file! ")
            print(" then run : gb_generator io_file\n ")

        except:
            print("Your input sigma is wrong!")
            sys.exit()

        return

if __name__ == '__main__':
    main()
