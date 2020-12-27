
# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
import numpy.linalg as la
from scipy.spatial import Delaunay
import matplotlib.tri as tri
from matplotlib.tri import TriRefiner, UniformTriRefiner
import matplotlib.pyplot as plt
import scipy
from scipy.special import binom as sc_binom

# ...
def barycentric_coords(vertices, point):
    T = (np.array(vertices[:-1])-vertices[-1]).T
    v = np.dot(la.inv(T), np.array(point)-vertices[-1])
    v.resize(len(vertices))
    v[-1] = 1-v.sum()
#    print v
    return v
# ...

class BezierTriangulationConstructor(object):
    def __init__(self, nodes, triangles=None, degree=None, control=None, weights=None, method="uniform"):
        # if triangles not given, then construct a Delaunay triangulation
        # creates a triangular bezier patch given the 3 summets of a triangle i and a total degree
        # A = vertices[i,0,:]
        # B = vertices[i,1,:]
        # C = vertices[i,2,:]
        self._array     = None
        self._local_triangles = []
        self._all_triangles = []
        self._n_patchs  = None
        self._degree = 1
        if degree is not None:
            self._degree = degree

        self._triangulation = tri.Triangulation(nodes[:,0], nodes[:,1], triangles)

        self._vertices = np.zeros((len(self._triangulation.triangles), 3, 2))
        self._vertices[:,:,0] = self._triangulation.x[self._triangulation.triangles]
        self._vertices[:,:,1] = self._triangulation.y[self._triangulation.triangles]

        if degree is not None:
            self._create_b_net(self._degree, self._vertices, control, weights, method)

        # strating from here self.degree, self.vertices, self.b_net must be initialized
        self._n_ptachs = self.vertices.shape[0]

        self._init_local_id()

        self._create_triangulation()
        self._update()

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
        return int((np.sqrt(8*self.shape+1) -1)//2) - 1

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
    def local_triangles(self):
        return self._local_triangles

    @property
    def all_triangles(self):
        return self._all_triangles

    @property
    def unique_points(self):
        return self._unique_points

    @property
    def unique_points_indices(self):
        return self._unique_points_idx

    @property
    def unique_triangles(self):
        return self._unique_triangles

    @property
    def triangulation(self):
        """
        returns the initial triangulation, without the B-nets
        """
        return self._triangulation

    @property
    def triangles_ancestors(self):
        return self._triangles_ancestors

    @property
    def edge(self, face):
        return 1

    def local_id(self, i, j, k):
        return self._local_id[i,j,k]

    def _init_local_id(self):
        degree = self.degree
        self._local_id = np.zeros((degree+1, degree+1, degree+1) \
                                   , dtype=np.int32)

        n = (degree+1)*(degree+2)//2
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
        n = (degree+1)*(degree+2)//2
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
        self._local_triangles = []
        self._all_triangles = []
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
            loc_triangles = np.array(loc_triangles)
            self._local_triangles.append(loc_triangles)
            n_tri = T_id * self.shape
            self._all_triangles.append(loc_triangles+n_tri)
        self._local_triangles = np.array(self._local_triangles)
        self._all_triangles = np.array(self._all_triangles)

    def _update(self):
        """
        computes unique vertices for the B-net
        """
        # ... points treatment
        points = self.points
        x_size = np.asarray(points.shape[:-1]).prod()
        A = np.zeros((x_size, points.shape[-1]))
        for i in range(0, points.shape[-1]):
            A[:, i] = points[...,i].reshape(x_size)
        a = np.ascontiguousarray(A)
        unique_a, idx = np.unique(a.view([('', a.dtype)]*a.shape[1]), return_index=True)
        self._unique_points = unique_a.view(a.dtype).reshape((unique_a.shape[0],a.shape[1]))
        self._unique_points_idx = idx

        # ... reverse indices
#        print idx
        reverse_idx = -np.ones(x_size, dtype=np.int32)
#        print ">>> Begin"
#        print a
        for i in range(0, x_size):
            itemindex = np.where(idx==i)[0]
#            print i, itemindex, a[i]
            if len(itemindex) > 0:
                reverse_idx[i] = itemindex[0]
            else:
                ll_condition = np.logical_and(a[i,0]== a[:, 0], a[i,1]== a[:, 1])
                k = np.where(ll_condition)
                reverse_idx[i] = reverse_idx[k[0][0]]
#        print "<<< End"
#        print reverse_idx

        # ... triangles treatment
        triangles = self.all_triangles
        ntri     = triangles.shape[0]
        ntri_loc = triangles.shape[1]
        nvert    = triangles.shape[2] # should be equal to 3
#        x_size = np.asarray(triangles.shape[:-1]).prod()
        x_size = ntri * ntri_loc
        B = np.zeros((x_size, nvert), dtype=np.int32)
        for i in range(0, triangles.shape[-1]):
            B[:, i] = reverse_idx[triangles[...,i].reshape(x_size)]
        self._unique_triangles = B

        triangles_ancestors = np.zeros((ntri, ntri_loc), dtype=np.int32)
        for T_id in range(0, ntri):
            triangles_ancestors[T_id,:] = T_id
        self._triangles_ancestors = triangles_ancestors.reshape(x_size)
#        print "ancestors"
#        print self._triangles_ancestors


class UniformBezierTriRefiner(TriRefiner):
    """
    Uniform Bezier mesh refinement

    Parameters
    ----------
    triangulation : :class:`~matplotlib.tri.Triangulation`
                     The encapsulated triangulation (to be refined)
    """
    def __init__(self, triangulation):
        TriRefiner.__init__(self, triangulation)

    def refine_triangulation(self, degree=2, ancestors=False):
        triang = self._triangulation
        x = triang.x ; y = triang.y
        triangles = triang.triangles
        nodes = np.zeros((x.shape[0],2))
        nodes [:,0] = x ; nodes [:,1] = y
        bzr = BezierTriangulationConstructor(nodes, triangles=triangles, degree=degree)

        triangles = bzr.unique_triangles
        points    = bzr.unique_points

        x = points[:,0] ; y = points[:,1]
        if not ancestors:
            return tri.Triangulation(x,y,triangles)
        else:
            return tri.Triangulation(x,y,triangles), bzr.triangles_ancestors



