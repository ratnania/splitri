# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from splitri.gallery.examples import collela, domain_1, domain_2, annulus
from splitri.bezier import Bezier

def test_1():
    L = 2.
    n = 30
    degree = 3

    triang = collela(n)
#    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    triangles = []
    for j in range(0,n-1):
        for i in range(0,n-1):
            I1 = i+j*n ; I2 = i+1+j*n ; I3 = i+1+(j+1)*n
            T = [I1,I2,I3]
            triangles.append(T)

            I1 = i+j*n ; I2 = i+(j+1)*n ; I3 = i+1+(j+1)*n
            T = [I1,I2,I3]
            triangles.append(T)

    Bzr = Bezier(degree, triang.x, triang.y, triangles)

    P = [0.5, 0.35]
    plt.plot(P[0], P[1], "or")
    T_id = Bzr.find_simplex(P)

    plt.plot(Bzr.triang.x[Bzr.triang.triangles[T_id]] \
             , Bzr.triang.y[Bzr.triang.triangles[T_id]], "og")

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()


def test_2():
    degree = 3

    triang = domain_1()

    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

#    P = [0.5, 0.35]
#    plt.plot(P[0], P[1], "or")
#    T_id = Bzr.find_simplex(P)
#
#    plt.plot(Bzr.triang.x[Bzr.triang.triangles[T_id]] \
#             , Bzr.triang.y[Bzr.triang.triangles[T_id]], "og")

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()


def test_3():
    degree = 4

    triang = domain_2()

    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

#    P = [0.5, 0.35]
#    plt.plot(P[0], P[1], "or")
#    T_id = Bzr.find_simplex(P)
#
#    plt.plot(Bzr.triang.x[Bzr.triang.triangles[T_id]] \
#             , Bzr.triang.y[Bzr.triang.triangles[T_id]], "og")

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()


def test_4():
    degree = 3
    min_radius = 0.25

    triang = annulus(n_radii = 8, n_angles = 36, min_radius = min_radius, max_radius = 0.95)

#    # Mask off unwanted triangles.
#    xmid = triang.x[triang.triangles].mean(axis=1)
#    ymid = triang.y[triang.triangles].mean(axis=1)
#
#    mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
#    triang.set_mask(mask)

    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    # Mask off unwanted triangles.
    triang = Bzr.triang
    xmid = triang.x[triang.triangles].mean(axis=1)
    ymid = triang.y[triang.triangles].mean(axis=1)

    mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triang.set_mask(mask)

    # Mask off unwanted triangles.
    triang = Bzr.triang_ref
    xmid = triang.x[triang.triangles].mean(axis=1)
    ymid = triang.y[triang.triangles].mean(axis=1)

    mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triang.set_mask(mask)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()

#######################################################
if __name__ == "__main__":
    test_1()
#    test_2()
#    test_3()
#    test_4()
