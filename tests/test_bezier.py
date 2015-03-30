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
    n = 3
    degree = 3

    triang = collela(n)
    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    x = Bzr.triang_ref.x
    y = Bzr.triang_ref.y

    Bzr._b_coeff = x**2+y**2

    npts = 100
    pts = np.random.random((npts,2))
    # pts <=1. but we need pts.sum() <=1.
    # one way to do it is the following
    pts[:,1] *= 1.-pts[:,0]

    positions = np.zeros((npts*Bzr.triang.triangles.shape[0],2))
    values    = np.zeros(npts*Bzr.triang.triangles.shape[0])
    i_pos = 0
    for T_id in range(0, Bzr.triang.triangles.shape[0]):
#        print ">>>>>>>>>>>>>> T ", T_id
        for i in range(0, npts):
            P = pts[i,:]
            if P.sum() < 1.:
#                print "pts", P
                v,position = Bzr.evaluate_on_triangle(P, T_id)
                positions[i_pos,:] = position
                values[i_pos]    = v
                i_pos += 1
#                plt.plot(position[0],position[1],".k")
#                print v, position

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")

    # create a triangulation for the computed positions
    triang_new = tri.Triangulation(positions[:,0],positions[:,1])
    plt.tripcolor(triang_new, values, shading='gouraud', cmap=plt.cm.rainbow)
    plt.triplot(triang_new, '-', lw=0.5, color="green")

#    xlim = [-1.2,1.2] ; ylim = [-1.2,1.2]
#    plt.xlim(*xlim)  ; plt.ylim(*ylim)
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
