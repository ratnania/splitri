# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
import numpy.linalg as la
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import scipy
from scipy.special import binom as sc_binom
from bezier_patch import *

def make_triangle_1(degree):
    vertices = np.zeros((3,2))
    #           A                B                     C
    vertices[0,0] = 0.0 ; vertices[1,0] = -1.0 ; vertices[2,0] = 1.0
    vertices[0,1] = 1.0 ; vertices[1,1] =  0.0 ; vertices[2,1] = 0.0
    return bezier_patch(degree, vertices)

def make_triangle_2():
    T = make_triangle_1(3)

    ids = T._ids_on_edges[0]
    points = T.array[ids,...]
    points[1:-1,1] -= 0.3
    T.array[ids,...] = points

    return T






#########################################################
if __name__ == "__main__":
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

#    for i in range(0, B_net.points.shape[0]):
#        plt.plot(B_net.points[i,0], B_net.points[i,1], "or")
    plt.triplot(B_net.points[:,0], B_net.points[:,1], B_net.triangles, 'b-')

#    for i in range(0, C_net.points.shape[0]):
#        plt.plot(C_net.points[i,0], C_net.points[i,1], "og")
#    plt.triplot(C_net.points[:,0], C_net.points[:,1], C_net.triangles, 'y-')

    crv_0 = A_net.extract_edge(0)
    crv_1 = A_net.extract_edge(1)
    crv_2 = A_net.extract_edge(2)

    t = np.linspace(0.,1.,10)
    for crv, col in zip([crv_0, crv_1, crv_2], ["r", "g", "y"]):
        P = crv(t)
        for i in range(0, P.shape[0]):
            plt.plot(P[i,0], P[i,1], "o"+col)

#    for i in range(0, P.shape[0]):
#        plt.plot(P[i,0], P[i,1], "or")

#    for i in range(0, P.shape[0]):
#        plt.plot(P[i,0], P[i,1], "og")

#    for i in range(0, P.shape[0]):
#        plt.plot(P[i,0], P[i,1], "oy")


#    plt.title('Bezier-net and its triangulation')
#    plt.savefig("bezier_triangle.png")

    plt.show()



