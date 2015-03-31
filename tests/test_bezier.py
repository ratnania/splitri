# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from splitri.gallery.examples import square,collela, domain_1, domain_2, annulus
from splitri.bezier import Bezier
from splitri.utils.triangle import barycentric_coords
from matplotlib.tri import Triangulation, UniformTriRefiner
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def test_1():
    L = 2.
    n = 10
    degree = 3

#    triang = collela(n)
    triang = square(n)

    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    x = Bzr.triang_ref.x
    y = Bzr.triang_ref.y

    Bzr._b_coeff = x**2+y**2

#    npts = 100
#    pts = np.random.random((npts,2))
#    pts = 2*pts-1.

#    old_settings = np.seterr(all='ignore')  #seterr to known value
#    np.seterr(over='raise')

    from splitri.boxsplines import triangulation_I
    n_I = 13 ; degree_I = 3 # degree
    tri_I = triangulation_I(n_I,degree_I,xmin=-1.,xmax=1.)

    triang_I = tri_I.triang
#    triang_I = tri_I.triang_ref
    npts = triang_I.x.shape[0]

    positions_x = []
    positions_y = []
    values    = []
    for i in range(0, npts):
        P = np.array([triang_I.x[i], triang_I.y[i]])

        # find in which triangle it belongs
        T_id = Bzr.find_simplex(P)
        if T_id is not None:
            # compute barycentric coordinates
            vertices = np.zeros((3,2))
            vertices[:,0] = Bzr.triang.x[Bzr.triang.triangles[T_id]]
            vertices[:,1] = Bzr.triang.y[Bzr.triang.triangles[T_id]]

            C = barycentric_coords(vertices, P)

            # evaluate the bezier surface on the point C within the triangle T
            v,position = Bzr.evaluate_on_triangle(C, T_id)
            positions_x.append(position[0])
            positions_y.append(position[1])
            values.append(v)

#    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
#    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")

    # create a triangulation for the computed positions
#    triang_new = tri.Triangulation(positions_x,positions_y)
    triang_new = triang_I
    plt.triplot(triang_new, '-', lw=0.5, color="green")

    refiner = UniformTriRefiner(triang_new)
    triang_new, values= refiner.refine_field(values, subdiv=3)

    values = np.array(values)
    print values.min(), values.max()

    plt.tripcolor(triang_new, values, shading='gouraud', cmap=plt.cm.rainbow)


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
