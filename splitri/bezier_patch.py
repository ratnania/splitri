# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
import numpy.linalg as la
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import scipy
from scipy.special import binom as sc_binom

from .bernstein import *

# ...
def barycentric_coords(vertices, point):
    T = (np.array(vertices[:-1])-vertices[-1]).T
    v = np.dot(la.inv(T), np.array(point)-vertices[-1])
    v.resize(len(vertices))
    v[-1] = 1-v.sum()
    return v
# ...

class bezier_patch(object):
    def __init__(self, degree=None, vertices=None, control=None, weights=None, method="uniform"):
        # creates a triangular bezier patch given the 3 summets of a triangle i and a total degree
        # A = vertices[i,0,:]
        # B = vertices[i,1,:]
        # C = vertices[i,2,:]
        self._array     = None
        self._triangles = []
        self._n_patchs  = None

        if vertices is not None:
            self._vertices = vertices
            if degree is not None:
                self._create_b_net(degree, vertices, control, weights, method)

        # strating from here self.degree, self.vertices, self.b_net must be initialized
        self._bernstein = bernstein(self.degree)
        self._n_ptachs = self.vertices.shape[0]

        self._init_local_id()

        self._create_triangulation()

    @property
    def dim(self):
        # total dim = dim + 1 (control) + 1 (weights/ not yet used)
        return self._array.shape[-1] - 1 - 1

    @property
    def vertices(self):
        return self._vertices

    @property
    def degree(self):
        # TODO must take the nearest integer rather than int
        return int((np.sqrt(8*self.shape+1) -1)/2) - 1

    @property
    def n_patchs(self):
        return self._array.shape[0]

    @property
    def shape(self):
        return self._array.shape[1]

    @property
    def array(self):
        return self._array

    @property
    def points(self):
        return self._array[...,:self.dim]

    @property
    def control(self):
        return self._array[:,:,self.dim]

    @property
    def weights(self):
        return self._array[:,:,self.dim+1]

    @property
    def triangles(self):
        return self._triangles

    @property
    def edge(self, face):
        return 1

    def local_id(self, i, j, k):
        return self._local_id[i,j,k]

    def _init_local_id(self):
        degree = self.degree
        self._local_id = np.zeros((degree+1, degree+1, degree+1) \
                                   , dtype=np.int32)

        n = (degree+1)*(degree+2)/2
        i_pos = 0
        for i in range(0, degree+1):
            for j in range(0, degree+1-i):
                k = degree - i - j
                self._local_id[i,j,k] = i_pos
                i_pos += 1

        self._ids_on_edges = []
        # face 0 : i = 0
        i = 0
        ids_on_edge = []
        for j in range(0, degree+1):
            k = degree - j
            ids_on_edge.append(self._local_id[i,j,k])
        self._ids_on_edges.append(ids_on_edge)

        # face 1 : j = 0
        j = 0
        ids_on_edge = []
        for i in range(0, degree+1):
            k = degree - i
            ids_on_edge.append(self._local_id[i,j,k])
        self._ids_on_edges.append(ids_on_edge)

        # face 2 : k = 0
        k = 0
        ids_on_edge = []
        for i in range(0, degree+1):
            j = degree - i
            ids_on_edge.append(self._local_id[i,j,k])
        self._ids_on_edges.append(ids_on_edge)

    def _create_b_net(self, degree, vertices, control, weights, method):
        n = (degree+1)*(degree+2)/2
        dim = vertices.shape[-1]
        n_patchs = vertices.shape[0]
        self._array = np.zeros((n_patchs,n,dim+1+1))
        if method == "uniform":
            for T_id in range(0, n_patchs):
                b_net = np.zeros((n,dim+1+1))
                # set weights to 1 by defualt
                b_net[:, dim+1] = 1.0
                i_pos = 0
                for i in range(0, degree+1):
                    for j in range(0, degree+1-i):
                        k = degree - i - j

                        b_net[i_pos, :dim] = i * vertices[T_id,0,:dim] \
                                           + j * vertices[T_id,1,:dim] \
                                           + k * vertices[T_id,2,:dim]

                        if control is not None:
                            b_net[i_pos, dim] = control[T_id,i_pos]

                        if weights is not None:
                            b_net[i_pos, dim+1] = weights[T_id,i_pos]

                        i_pos += 1

                b_net /= degree
                self._array[T_id, ...] = b_net

    def _create_triangulation(self):
        self._triangles = []
        for T_id in range(0, self.n_patchs):
            loc_triangles = []
            # top triangles
            for i in range(0, self.degree+1-1):
                for j in range(0, self.degree+1-i-1):
                    k = self.degree - i - j
                    I1 = self.local_id(i,j,k)
                    I2 = self.local_id(i,j+1,k-1)
                    I3 = self.local_id(i+1,j,k-1)
    #                print "UP   ", [I1, I2, I3]
                    loc_triangles.append([I1, I2, I3])
            # bottom triangles
            for i in range(1, self.degree+1):
                for j in range(0, self.degree+1-i-1):
                    k = self.degree - i - j
                    I1 = self.local_id(i,j,k)
                    I2 = self.local_id(i,j+1,k-1)
                    I3 = self.local_id(i-1,j+1,k)
    #                print "DOWN ", [I1, I2, I3]
                    loc_triangles.append([I1, I2, I3])
            self._triangles.append(np.array(loc_triangles))
        self._triangles = np.array(self._triangles)

    def set_b_coefficients(self, control):
        self._array[...,self.dim] = control

    def copy(self):
        bzr = bezier_patch.__new__(type(self))
        bzr._array = self.array.copy()
        bzr._triangles = self.triangles.copy()
        return bzr

    def clone(self):
        bzr = bezier_patch.__new__(type(self))
        bzr._array = self.array
        bzr._triangles = self.triangles
        return bzr

    def find_simplex(self, xyz):
        """
        returns the id of the triangle that contains the point xyz
        """
#        c = barycentric_coords(self.vertices, xyz)
        # TODO must be done using numpy arrays without a loop
        for i in range(0, self.vertices.shape[0]):
            c = barycentric_coords(self.vertices[i,:], xyz)
            if np.array(c).prod() >= 0:
                return i
        return None

    def extract_edge(self, T_id, edge_id):
        ids = self._ids_on_edges[edge_id]
        points = self.array[T_id, ids,:self.dim]
        from igakit.nurbs import NURBS
        degree= self.degree
        U = [0.]*(degree+1) + [1.]*(degree+1)
        nrb = NURBS([U], points)
        return nrb

    def translate(self, displ, axis=None):
        displ = np.asarray(displ, dtype='d')
        assert displ.ndim in (0, 1)
        if displ.ndim > 0 and axis is None:
            assert displ.size <= 2
            t = displ
        else:
            t = np.zeros(2, dtype='d')
            if axis is None:
                t[0] = displ
            else:
                t[axis] = displ

#        self._array += displ
        self._array[:,:,:self.dim] += displ
        return self

    def scale(self, scale, axis=None):
        scale = np.asarray(scale, dtype='d')
        assert scale.ndim in (0, 1)
        if scale.ndim > 0 and axis is None:
            assert scale.size <= 2
            t = scale
        else:
            t = np.zeros(2, dtype='d')
            if axis is None:
                t[0] = scale
            else:
                t[axis] = scale
#        self._array *= scale
        self._array[:,:,:self.dim] *= scale
        return self

#    def move(self, displ, axis=None):
#        t = transform().move(displ, axis)
#        return self.transform(t)
#
#    def rotate(self, angle, axis=2):
#        t = transform().rotate(angle, axis)
#        return self.transform(t)

    #

    def subdivise(self):
        pass

    def elevate(self, times):
        pass

    def __call__(self, xy):
        """
        evaluates the surface using bernstein polynomial
        """
        A = np.array(xy)
        T_id = self.find_simplex(A)
        vertices = self.vertices[T_id,...]
        coeffs = self.control[T_id,...]
        c = barycentric_coords(vertices, A)

        value = 0.
        i_pos = 0
        for i in range(0, self.degree+1):
            for j in range(0, self.degree+1-i):
                k = self.degree - i - j
                bern = self._bernstein([i,j,k], c)
                value += bern * coeffs[i_pos]
                i_pos += 1
        return value

    def plot(self):
        for T_id in range(0, self.n_patchs):
            points = self.points[T_id,...]
    #        for i in range(0, points.shape[0]):
    #            plt.plot(points[i,0], points[i,1], "or")
            triangles = self.triangles[T_id, ...]
            plt.triplot(points[:,0], points[:,1], triangles, 'b-')

            z = []
            for i in range(0, points.shape[0]):
                z.append(self(points[i,:]))
            z = np.array(z)
            plt.tricontourf(points[:,0], points[:,1], triangles, z)
