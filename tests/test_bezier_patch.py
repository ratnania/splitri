# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
import numpy.linalg as la
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import scipy
from scipy.special import binom as sc_binom
from splitri.bezier_patch import *
from splitri.bernstein import *

def make_triangle_1(degree):
    vertices = np.zeros((1,3,2))
    #           A                B                     C
    vertices[0,0,0] = 0.0 ; vertices[0,1,0] = -1.0 ; vertices[0,2,0] = 1.0
    vertices[0,0,1] = 1.0 ; vertices[0,1,1] =  0.0 ; vertices[0,2,1] = 0.0
    return bezier_patch(degree, vertices)

def make_triangle_2(degree):
    T = make_triangle_1(degree)

    ids = T._ids_on_edges[0]
    points = T.array[0,ids,...]
    points[1:-1,1] -= 0.3
    T.array[0,ids,...] = points

    return T

def make_triangles_1(degree):
    vertices = np.zeros((2,3,2))
    #           A                B                     C
    # ... Triangle id = 0
    T_id = 0
    vertices[T_id,0,0] = 0.0 ; vertices[T_id,1,0] = -1.0 ; vertices[T_id,2,0] = 1.0
    vertices[T_id,0,1] = 1.0 ; vertices[T_id,1,1] =  0.0 ; vertices[T_id,2,1] = 0.0
    # ... Triangle id = 1
    T_id = 1
    vertices[T_id,0,0] = 1.0 ; vertices[T_id,1,0] = -1.0 ; vertices[T_id,2,0] = 0.0
    vertices[T_id,0,1] = 0.0 ; vertices[T_id,1,1] =  0.0 ; vertices[T_id,2,1] = -1.0
    return bezier_patch(degree, vertices)

def make_triangles_square(degree):
    # ...
    n = 10
    u = np.linspace(0.,1.,n)
    v = np.linspace(0.,1.,n)

    U,V = np.meshgrid(u,v)

    u = U.reshape(U.size)
    v = V.reshape(V.size)
    # ...

    # ...
    points = np.zeros((U.size, 2))
    points[:,0] = u
    points[:,1] = v
    # ...

    # ... create the tesselation, triangulation
    tri = Delaunay(points)
    n_triangles = tri.simplices.shape[0]
    vertices = np.zeros((n_triangles, 3, 2))
    for T_id in range(0, n_triangles):
        ids = tri.simplices[T_id,:]
        vertices[T_id, :, :] = points[ids, :]
    # ...

    return bezier_patch(degree, vertices)

def make_triangles_2(degree):
    # ...
    n = 10
    u = np.linspace(0.,1.,n)
    v = np.linspace(0.,1.,n)

    U,V = np.meshgrid(u,v)

    u = U.reshape(U.size)
    v = V.reshape(V.size)
    # ...

    # ...
    eps = 0.1
    k1 = 1. ; k2 = 1.
    X = 2*(u + eps * sin(2*pi*k1*u) * sin(2*pi*k2*v)) -1.
    Y = 2*(v + eps * sin(2*pi*k1*u) * sin(2*pi*k2*v)) - 1.
    # ...

    # ...
    points = np.zeros((X.size, 2))
    points[:,0] = X
    points[:,1] = Y
    # ...

    # ... create the tesselation, triangulation
    tri = Delaunay(points)
    n_triangles = tri.simplices.shape[0]
    vertices = np.zeros((n_triangles, 3, 2))
    for T_id in range(0, n_triangles):
        ids = tri.simplices[T_id,:]
        vertices[T_id, :, :] = points[ids, :]
    # ...

    return bezier_patch(degree, vertices)

def make_triangles_3(degree):
    points = np.zeros((4,2))
    points[0,:] = [0.0, 1.0] # A
    points[1,:] = [-1., 0.0] # B
    points[2,:] = [1.0, 0.0] # C
    points[3,:] = [0.0, -1.] # D

    # ... create the tesselation, triangulation
    tri = Delaunay(points)
    n_triangles = tri.simplices.shape[0]
    vertices = np.zeros((n_triangles, 3, 2))
    for T_id in range(0, n_triangles):
        ids = tri.simplices[T_id,:]
        vertices[T_id, :, :] = points[ids, :]
    # ...

    return bezier_patch(degree, vertices)

def test_1():
#    A_net = make_triangle_1(3)
    A_net = make_triangle_2()
    B_net = A_net.copy()

#    C_net = A_net.clone()
#
#    C_net.translate([1.0,1.0])
#    C_net.scale([0.5,0.5])
#
#    print C_net.degree

    xlim = [-2.2, 2.2]
    ylim = [-2.2, 2.2]

    plt.figure()
#    plt.axis('off')
#    plt.gca().set_aspect('equal')
    plt.xlim(*xlim)
    plt.ylim(*ylim)

    for T_id in range(0, A_net.n_patchs):
        points = A_net.points[T_id,...]
#        for i in range(0, points.shape[0]):
#            plt.plot(points[i,0], points[i,1], "or")
        triangles = A_net.triangles[T_id, ...]
        plt.triplot(points[:,0], points[:,1], triangles, 'b-')

#    for i in range(0, B_net.points.shape[0]):
#        plt.plot(B_net.points[i,0], B_net.points[i,1], "or")
#    plt.triplot(B_net.points[:,0], B_net.points[:,1], B_net.triangles, 'b-')

#    for i in range(0, C_net.points.shape[0]):
#        plt.plot(C_net.points[i,0], C_net.points[i,1], "og")
#    plt.triplot(C_net.points[:,0], C_net.points[:,1], C_net.triangles, 'y-')

    T_id = 0
    crv_0 = A_net.extract_edge(T_id, 0)
    crv_1 = A_net.extract_edge(T_id, 1)
    crv_2 = A_net.extract_edge(T_id, 2)

    t = np.linspace(0.,1.,50)
    for crv, col in zip([crv_0, crv_1, crv_2], ["r", "g", "y"]):
        P = crv(t)
        for i in range(0, P.shape[0]):
            plt.plot(P[i,0], P[i,1], ".-"+col)

#    for i in range(0, P.shape[0]):
#        plt.plot(P[i,0], P[i,1], "or")

#    for i in range(0, P.shape[0]):
#        plt.plot(P[i,0], P[i,1], "og")

#    for i in range(0, P.shape[0]):
#        plt.plot(P[i,0], P[i,1], "oy")


#    plt.title('Bezier-net and its triangulation')
#    plt.savefig("bezier_triangle.png")

    plt.show()


def test_2():
    bzr = make_triangles_1(3)
    bzr_1 = bzr.copy()

    print bzr.find_simplex([0.5,0.5])

    xlim = [-2.2, 2.2]
    ylim = [-2.2, 2.2]

    plt.figure()
#    plt.axis('off')
#    plt.gca().set_aspect('equal')
    plt.xlim(*xlim)
    plt.ylim(*ylim)

    for T_id in range(0, bzr.n_patchs):
        points = bzr.points[T_id,...]
#        for i in range(0, points.shape[0]):
#            plt.plot(points[i,0], points[i,1], "or")
        triangles = bzr.triangles[T_id, ...]
        plt.triplot(points[:,0], points[:,1], triangles, 'b-')

    plt.show()

def test_3():
#    bzr = make_triangles_2(5)
    bzr = make_triangles_square(5)
    bzr_1 = bzr.copy()

    xlim = [-.2, 1.2]
    ylim = [-.2, 1.2]

    plt.figure()
#    plt.axis('off')
#    plt.gca().set_aspect('equal')
    plt.xlim(*xlim)
    plt.ylim(*ylim)

    for T_id in range(0, bzr.n_patchs):
        points = bzr.points[T_id,...]
#        for i in range(0, points.shape[0]):
#            plt.plot(points[i,0], points[i,1], "or")
        triangles = bzr.triangles[T_id, ...]
        plt.triplot(points[:,0], points[:,1], triangles, 'b-')

    plt.show()

def test_4():
    degree = 5
#    bzr = make_triangles_3(degree)
    bzr = make_triangles_square(degree)
#    b_coeff = np.random.rand(bzr.shape)
    x = bzr.points[...,0]
    y = bzr.points[...,1]
    sin = np.sin ; pi = np.pi
#    b_coeff = np.ones_like(bzr.points[...,0])
#    b_coeff = sin(2*pi*x)*sin(2*pi*y)
    b_coeff = x**2 + y**2

    bzr.set_b_coefficients(b_coeff)


    plt.figure()
#    plt.axis('off')
#    plt.gca().set_aspect('equal')
#    xlim = [-2.2, 2.2]
#    ylim = [-2.2, 2.2]
#    plt.xlim(*xlim)
#    plt.ylim(*ylim)

#    bzr.plot(nlevel=20, vmin=-0.1, vmax=2.1)
    bzr.plot()
#    plt.colorbar()
    plt.show()

#########################################################
if __name__ == "__main__":
#    test_1()
#    test_2()
#    test_3()
    test_4()
