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

def make_triangles_6(degree):

    x =np.array([-1.1288,-0.27786,0.80753,1.0593,-0.1563,-0.62518,-0.95861,-0.78842,-0.61823,-0.44805,-0.096961,0.083936,0.26483,0.44573,0.62663,0.85789,0.90825,0.95861,1.009,0.85673,0.65412,0.45152,0.24891,0.04631,-0.27352,-0.39074,-0.50796,-0.79305,-0.96093,0.093606,-0.70378,0.72463,-0.27503,0.64406,-0.30976,0.40348,0.28319,-0.10986,-0.073193,0.87604,-0.88885,0.19124,-0.00036351,-0.51538,-0.3409,0.68238,0.43689,-0.6176,0.54328,-0.079635,0.31319,0.73076,-0.79277,0.87668,-0.20567,-0.21595,0.11589,0.26013,0.32212,0.54986,0.45791,0.12746,-0.44664,-0.28559,0.11883,0.061646,-0.50891,-0.48716,-0.62684,0.57669,0.74722,0.81603,0.37258,0.22964,-0.41324,-0.1382,-0.37681,-0.035599,0.037716,-0.068816,-0.22796,-0.060578,-0.43952,-0.20434])
    y =np.array([0.11288,0.68162,0.23444,-0.60781,-0.75543,-0.29088,0.22663,0.34038,0.45412,0.56787,0.60709,0.53256,0.45803,0.3835,0.30897,0.065991,-0.10246,-0.27091,-0.43936,-0.63242,-0.65702,-0.68162,-0.70622,-0.73082,-0.63929,-0.52315,-0.40702,-0.1563,-0.021708,-0.11758,0.14118,-0.37025,0.45932,0.091961,0.11512,-0.16654,0.13428,-0.36803,0.3966,-0.48949,0.13423,-0.40068,0.1352,0.31481,-0.20473,-0.21478,0.01804,-0.055294,-0.48544,-0.56999,0.29215,-0.52686,0.0078785,-0.36062,0.26627,-0.065918,0.28055,-0.050238,-0.53119,-0.28196,0.20482,-0.56317,0.41544,-0.35988,0.061395,-0.29014,0.14657,-0.18565,0.27854,-0.10593,-0.083011,-0.23355,-0.34932,-0.22943,-0.043161,0.11161,0.2849,-0.010632,-0.43886,-0.18259,-0.49244,0.23716,-0.32913,-0.23735])

    t1 =np.array([7,28,8,9,11,10,2,12,14,16,15,3,17,18,20,19,4,21,22,34,25,23,5,26,13,29,31,33,47,41,44,40,24,68,48,27,58,8,66,50,49,11,54,60,55,39,55,49,46,1,48,40,64,58,51,57,56,56,64,16,47,57,47,67,6,60,73,66,59,12,51,20,32,31,28,32,71,65,63,76,68,76,37,78,36,59,22,32,66,37,14,62,23,9,35,80,50,37,30,36,38,64,31,67,45,67,31,34,36,70,34,32,17,42,49,30,42,35,48,39,35,33,44,30,43,50,42,38,30,25,38,43,55,26,45,45,38])
    t2 =np.array([1,6,7,8,2,9,10,11,13,3,14,15,16,17,4,18,19,20,21,15,5,22,24,25,12,28,8,10,34,29,9,19,23,6,6,26,30,31,30,24,32,33,18,36,35,33,33,21,32,29,31,32,38,37,13,39,35,45,26,34,37,43,36,35,27,46,36,38,42,39,37,40,49,41,48,40,46,43,44,56,45,43,51,56,47,49,49,46,42,47,51,42,59,44,55,56,38,57,58,58,50,45,48,48,56,44,67,47,60,46,70,54,71,59,60,66,73,67,68,55,56,63,67,65,76,62,66,66,78,50,64,57,76,64,68,64,80])
    t3 =np.array([41,48,41,69,33,63,33,39,51,34,61,34,71,72,40,54,40,52,49,61,50,59,50,81,57,53,41,63,61,53,69,54,62,83,68,83,74,69,80,62,60,39,72,73,76,55,77,52,72,41,53,52,84,65,57,82,75,84,81,71,58,65,70,77,83,70,74,79,62,57,61,52,52,53,53,54,72,78,77,78,75,82,57,80,58,73,59,60,74,61,61,79,62,63,77,84,81,65,65,74,79,83,67,75,75,69,69,70,70,71,71,72,72,73,73,74,74,75,75,82,76,77,77,78,78,79,79,80,80,81,81,82,82,83,83,84,84])

    tri = np.vstack((t1-1,t2-1,t3-1)).transpose()

    nodes = np.zeros((x.shape[0], 2))
    nodes[:,0] = x ; nodes[:,1] = y

    return bezier_patch(nodes, triangles=tri, degree=degree)

def make_triangles_square_boxspline_I(n,degree):
    t = np.linspace(0.,1.,n)
    X,Y=np.meshgrid(t,t)
    x = X.reshape(X.size)
    y = Y.reshape(Y.size)

    i_t = list(range(0,n))
    i_X,i_Y = np.array(np.meshgrid(i_t,i_t), dtype=np.int32)

    i_x = i_X.reshape(i_X.size)
    i_y = i_Y.reshape(i_Y.size)

    ids = i_x + i_y*n

    triangles = []
    for j in range(0,n-1):
        for i in range(0,n-1):
            I1 = i+j*n ; I2 = i+1+j*n ; I3 = i+1+(j+1)*n
            T = [I1,I2,I3]
            triangles.append(T)

            I1 = i+j*n ; I2 = i+(j+1)*n ; I3 = i+1+(j+1)*n
            T = [I1,I2,I3]
            triangles.append(T)

    nodes = np.zeros((x.shape[0], 2))
    nodes[:,0] = x ; nodes[:,1] = y

    triang = tri.Triangulation(x,y,triangles)
#    plt.triplot(x,y,triang.triangles)
#    xlim = [-.2, 1.2]
#    ylim = [-.2, 1.2]
#    plt.xlim(*xlim)
#    plt.ylim(*ylim)
#    plt.show()
#    import sys; sys.exit(0)

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

    print(bzr.find_simplex([0.5,0.5]))

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
    degree = 2
#    bzr = make_triangles_1(degree)
#    bzr = make_triangles_4(degree)
#    bzr = make_triangles_6(degree)
###    bzr = make_triangles_5(degree)
#    bzr = make_triangles_collela(20,degree)
#    bzr = make_triangles_square(3,degree)
    bzr = make_triangles_square_boxspline_I(6,degree)
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
    bzr.plot(show_triangles=True, show_values=True, show_triangulation=True)

    plt.colorbar()
    plt.show()

def test_boxspline_1():
    n = 4
    degree = 5
    bzr = make_triangles_square_boxspline_I(n,degree)
    bzr.update()

    x = bzr.points[...,0]
    y = bzr.points[...,1]

#    print bzr.find_simplex([0.75, 0.7])
#    print bzr.find_simplex([0.65, 0.5])

    b_coeff = np.zeros_like(x)
    print(b_coeff.shape)

#    b_coeff[:10,0:6] = [1.,1.,1.,1.,1.,1.]
#    b_coeff[0,0:6] = [1.,1.,1.,1.,1.,1.]

#    b_coeff[0:10,0] = 1.
#    plt.plot(x[0,0:6], y[0,0:6],"ob")

    b_coeff[10:20,0] = 1.
    plt.plot(x[10:20,:], y[10:20,:],"ob")

#    b_coeff[34,0:6] = [1.,1.,1.,1.,1.,1.]
#    b_coeff[34,0:6] = [1.,2.,1.,1.,1.,0.]
#    b_coeff[27,0:6] = [1.,2.,1.,1.,1.,0.]
#    b_coeff[27,0:6] = [0.,1.,1.,1.,2.,1.]

    bzr.set_b_coefficients(b_coeff)

    print(bzr.unique_points.shape)


#    plt.figure()
#    plt.axis('off')
#    plt.gca().set_aspect('equal')
#    xlim = [-2.2, 2.2]
#    ylim = [-2.2, 2.2]
#    plt.xlim(*xlim)
#    plt.ylim(*ylim)

#    bzr.plot(nlevel=20, vmin=-0.1, vmax=2.1)
    bzr.plot(show_triangles=True, show_values=True,
             show_triangulation=True,subdiv=3)

    plt.colorbar()
    plt.show()

#########################################################
if __name__ == "__main__":
#    test_1()
#    test_2()
#    test_3()
#    test_4()
    test_boxspline_1()

