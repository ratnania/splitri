# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
import numpy.linalg as la
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import scipy
from scipy.special import binom as sc_binom

# ...
targets = np.array([[.1,0.], [.9,.49], [.1,.6], [.4,.9]])
# ...

# ... pre-compute all binomial coef needed
degree = 4
allBinom = []
for d in range(0,degree+1):
    values = np.zeros(d+1)
    for i in range(0, d+1):
        values[i] = sc_binom(d,i)
    allBinom.append(values)

#for d in range(0, degree+1):
#    print(len(allBinom[d]))

def binom(n,i):
    return allBinom[n][i]
# ...

# ...
def barycentric_coords(vertices, point):
    T = (np.array(vertices[:-1])-vertices[-1]).T
    v = np.dot(la.inv(T), np.array(point)-vertices[-1])
    v.resize(len(vertices))
    v[-1] = 1-v.sum()
    return v
# ...

# ...
def bezier(x,y,i,j,n):
    A = np.array([x,y])
    t_id = tri.find_simplex(targets)[0]
    vertices = tri.points[tri.vertices[t_id,:]]
    c = barycentric_coords(vertices, A)
    t0 = c[0] ; t1 = c[1] ; t2 = c[2]
    k = n - j - i
    v = 1.
#    v = t0**i * t1**j * t2**k
    v *= binom(n,i)
    v *= binom(n-i,j)
    return v
# ...

class bezier_patch(object):
    def __init__(self, degree=None, vertices=None, control=None, method="uniform"):
        # creates a triangular bezier patch given the 3 summets of a triangle i and a total degree
        # A = vertices[i,0,:]
        # B = vertices[i,1,:]
        # C = vertices[i,2,:]
        self._array     = None
        self._triangles = []
        self._n_patchs  = None

        if control is None:
            assert degree is not None
            assert vertices is not None
#            assert vertices is not None

            self._vertices = vertices
            self._create_b_net(degree, vertices, method)

        # strating from here self.degree, self.vertices, self.b_net must be initialized
        self._n_ptachs = self.vertices.shape[0]

        self._init_local_id()

        self._create_triangulation()

    @property
    def dim(self):
        return self._array.shape[-1]

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
        return self._array

    @property
    def control(self):
        return self._array

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

    def _create_b_net(self, degree, vertices, method):
        n = (degree+1)*(degree+2)/2
        dim = vertices.shape[-1]
        n_patchs = vertices.shape[0]
        b_net = np.zeros((n,dim))
        self._array = np.zeros((n_patchs,n,dim))
        if method == "uniform":
            for T_id in range(0, n_patchs):
                i_pos = 0
                for i in range(0, degree+1):
                    for j in range(0, degree+1-i):
                        k = degree - i - j

                        b_net[i_pos, :] = i * vertices[T_id,0,:] \
                                        + j * vertices[T_id,1,:] \
                                        + k * vertices[T_id,2,:]

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

    def extract_edge(self, T_id, edge_id):
        ids = self._ids_on_edges[edge_id]
        points = self.array[T_id, ids,...]
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
        self._array += displ
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
        self._array *= scale
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

