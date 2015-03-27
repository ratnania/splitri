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
import math

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
    points = np.zeros((4,2))
    points[0,:] = [0.0, 1.0] # A
    points[1,:] = [-1., 0.0] # B
    points[2,:] = [1.0, 0.0] # C
    points[3,:] = [0.0, -1.] # D

    return bezier_patch(points, degree=degree)

def make_triangles_square(n,degree):
    # ...
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

    return bezier_patch(points, degree=degree)

def make_triangles_collela(n,degree):
    # ...
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

    return bezier_patch(points, degree=degree)

def make_triangles_4(degree):

    xy = np.asarray([
        [-0.101, 0.872], [-0.080, 0.883], [-0.069, 0.888], [-0.054, 0.890],
        [-0.045, 0.897], [-0.057, 0.895], [-0.073, 0.900], [-0.087, 0.898],
        [-0.090, 0.904], [-0.069, 0.907], [-0.069, 0.921], [-0.080, 0.919],
        [-0.073, 0.928], [-0.052, 0.930], [-0.048, 0.942], [-0.062, 0.949],
        [-0.054, 0.958], [-0.069, 0.954], [-0.087, 0.952], [-0.087, 0.959],
        [-0.080, 0.966], [-0.085, 0.973], [-0.087, 0.965], [-0.097, 0.965],
        [-0.097, 0.975], [-0.092, 0.984], [-0.101, 0.980], [-0.108, 0.980],
        [-0.104, 0.987], [-0.102, 0.993], [-0.115, 1.001], [-0.099, 0.996],
        [-0.101, 1.007], [-0.090, 1.010], [-0.087, 1.021], [-0.069, 1.021],
        [-0.052, 1.022], [-0.052, 1.017], [-0.069, 1.010], [-0.064, 1.005],
        [-0.048, 1.005], [-0.031, 1.005], [-0.031, 0.996], [-0.040, 0.987],
        [-0.045, 0.980], [-0.052, 0.975], [-0.040, 0.973], [-0.026, 0.968],
        [-0.020, 0.954], [-0.006, 0.947], [ 0.003, 0.935], [ 0.006, 0.926],
        [ 0.005, 0.921], [ 0.022, 0.923], [ 0.033, 0.912], [ 0.029, 0.905],
        [ 0.017, 0.900], [ 0.012, 0.895], [ 0.027, 0.893], [ 0.019, 0.886],
        [ 0.001, 0.883], [-0.012, 0.884], [-0.029, 0.883], [-0.038, 0.879],
        [-0.057, 0.881], [-0.062, 0.876], [-0.078, 0.876], [-0.087, 0.872],
        [-0.030, 0.907], [-0.007, 0.905], [-0.057, 0.916], [-0.025, 0.933],
        [-0.077, 0.990], [-0.059, 0.993]])
    x = xy[:, 0]*180/3.14159
    y = xy[:, 1]*180/3.14159

    triangles = np.asarray([
        [67, 66,  1], [65,  2, 66], [ 1, 66,  2], [64,  2, 65], [63,  3, 64],
        [60, 59, 57], [ 2, 64,  3], [ 3, 63,  4], [ 0, 67,  1], [62,  4, 63],
        [57, 59, 56], [59, 58, 56], [61, 60, 69], [57, 69, 60], [ 4, 62, 68],
        [ 6,  5,  9], [61, 68, 62], [69, 68, 61], [ 9,  5, 70], [ 6,  8,  7],
        [ 4, 70,  5], [ 8,  6,  9], [56, 69, 57], [69, 56, 52], [70, 10,  9],
        [54, 53, 55], [56, 55, 53], [68, 70,  4], [52, 56, 53], [11, 10, 12],
        [69, 71, 68], [68, 13, 70], [10, 70, 13], [51, 50, 52], [13, 68, 71],
        [52, 71, 69], [12, 10, 13], [71, 52, 50], [71, 14, 13], [50, 49, 71],
        [49, 48, 71], [14, 16, 15], [14, 71, 48], [17, 19, 18], [17, 20, 19],
        [48, 16, 14], [48, 47, 16], [47, 46, 16], [16, 46, 45], [23, 22, 24],
        [21, 24, 22], [17, 16, 45], [20, 17, 45], [21, 25, 24], [27, 26, 28],
        [20, 72, 21], [25, 21, 72], [45, 72, 20], [25, 28, 26], [44, 73, 45],
        [72, 45, 73], [28, 25, 29], [29, 25, 31], [43, 73, 44], [73, 43, 40],
        [72, 73, 39], [72, 31, 25], [42, 40, 43], [31, 30, 29], [39, 73, 40],
        [42, 41, 40], [72, 33, 31], [32, 31, 33], [39, 38, 72], [33, 72, 38],
        [33, 38, 34], [37, 35, 38], [34, 38, 35], [35, 37, 36]])

    return bezier_patch(xy, triangles=triangles, degree=degree)

def make_triangles_5(degree):
    # First create the x and y coordinates of the points.
    n_angles = 36
    n_radii = 8
    min_radius = 0.25
    radii = np.linspace(min_radius, 0.95, n_radii)

    angles = np.linspace(0, 2*math.pi, n_angles, endpoint=False)
    angles = np.repeat(angles[..., np.newaxis], n_radii, axis=1)
    angles[:, 1::2] += math.pi/n_angles

    x = (radii*np.cos(angles)).flatten()
    y = (radii*np.sin(angles)).flatten()
    z = (np.cos(radii)*np.cos(angles*3.0)).flatten()

    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triang = tri.Triangulation(x, y)

#    # Mask off unwanted triangles.
#    xmid = x[triang.triangles].mean(axis=1)
#    ymid = y[triang.triangles].mean(axis=1)
#    mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
#    triang.set_mask(mask)

    nodes = np.zeros((x.shape[0], 2))
    nodes[:,0] = x ; nodes[:,1] = y

    return bezier_patch(nodes, degree=degree)


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
#    bzr = make_triangles_collela(10,5)
    bzr = make_triangles_square(10,5)
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
#    bzr = make_triangles_1(degree)
#    bzr = make_triangles_4(degree)
###    bzr = make_triangles_5(degree)
    bzr = make_triangles_collela(20,degree)
#    bzr = make_triangles_square(10,degree)
#    b_coeff = np.random.rand(bzr.shape)
    x = bzr.points[...,0]
    y = bzr.points[...,1]
    sin = np.sin ; pi = np.pi
#    b_coeff = np.ones_like(bzr.points[...,0])
#    b_coeff = sin(2*pi*x) #*sin(2*pi*y)
    b_coeff = x**2 + y**2
#    b_coeff = 16.*x*(1-x)*y*(1-y)
#    b_coeff = np.zeros_like(x)
#    b_coeff[0,0:3,...] = 1.

    bzr.set_b_coefficients(b_coeff)
    bzr.update()

    plt.figure()
#    plt.axis('off')
#    plt.gca().set_aspect('equal')
#    xlim = [-2.2, 2.2]
#    ylim = [-2.2, 2.2]
#    plt.xlim(*xlim)
#    plt.ylim(*ylim)

#    bzr.plot(nlevel=20, vmin=-0.1, vmax=2.1)
    bzr.plot(show_triangles=True, show_values=False)
#    plt.colorbar()

#    triangles = bzr.unique_triangles
#    points = bzr.unique_points
#    plt.triplot(points[:,0], points[:,1], triangles, 'b-')

    plt.show()

#########################################################
if __name__ == "__main__":
#    test_1()
#    test_2()
#    test_3()
    test_4()

#def unique(a):
#    order = np.lexsort(a.T)
#    a = a[order]
#    diff = np.diff(a, axis=0)
#    ui = np.ones(len(a), 'bool')
#    ui[1:] = (diff != 0).any(axis=1)
#    return diff, ui, a[ui]
#
#def unique_rows(a):
#    a = np.ascontiguousarray(a)
#    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
#    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1])), unique_a
