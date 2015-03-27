# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
import numpy.linalg as la
from scipy.spatial import Delaunay
import matplotlib.tri as tri
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
        self._local_triangles = []
        self._all_triangles = []
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

    def update(self):
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
        print ">>> Begin"
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
        print "<<< End"
#        print reverse_idx

#        # ... reverse indices
#        reverse_idx = -np.ones(x_size, dtype=np.int32)
#        for i in range(0, x_size):
#            try:
#                reverse_idx[i] = list(idx).index(i)
#            except:
#                for j in range(0, a.shape[0]):
#                    if np.linalg.norm(a[i]-a[j]) < 1.e-7:
#                        reverse_idx[i] = reverse_idx[j]
#        print reverse_idx

        # ... triangles treatment
        triangles = self.all_triangles
        x_size = np.asarray(triangles.shape[:-1]).prod()
        B = np.zeros((x_size, triangles.shape[-1]), dtype=np.int32)
        for i in range(0, triangles.shape[-1]):
            B[:, i] = reverse_idx[triangles[...,i].reshape(x_size)]
        self._unique_triangles = B

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

    def plot(self, show_triangles=True, show_values=False):
        triangles = self.unique_triangles
        points = self.unique_points

        x = points[:,0]
        y = points[:,1]

        # Create the Triangulation; no triangles so Delaunay triangulation created.
        triang = tri.Triangulation(x, y, triangles)

        if show_values:
            z = np.zeros_like(x)
            for i in range(0, z.shape[0]):
                z[i] = self.__call__([x[i], y[i]])
            plt.tripcolor(triang, z, shading='gouraud', cmap=plt.cm.rainbow)

        if show_triangles:
            plt.triplot(points[:,0], points[:,1], triangles, lw=0.5)

#        # Mask off unwanted triangles.
#        xmid = x[triangles].mean(axis=1)
#        ymid = y[triangles].mean(axis=1)
#
#        zfaces = np.zeros_like(xmid)
#        for i in range(0, zfaces.shape[0]):
#            zfaces[i] = self.__call__([xmid[i], ymid[i]])
#        plt.tripcolor(x, y, triangles, facecolors=zfaces, shading='gouraud',
#                      cmap=plt.cm.rainbow) #, edgecolors='k')

#    def plot(self):
#        triangles = []
#        points    = []
#        for T_id in range(0, self.n_patchs):
#            for loc_T_id in range(0, self.all_triangles.shape[1]):
#                triangles.append(list(self.all_triangles[T_id, loc_T_id, ...]))
#            for loc_id in range(0, self.points.shape[1]):
#                points.append(list(self.points[T_id, loc_id, ...]))
#        points    = np.array(points)
#        triangles = np.array(triangles, dtype=np.int32)
#
#        x = points[:,0]
#        y = points[:,1]
#        xmid = x[triangles].mean(axis=1)
#        ymid = y[triangles].mean(axis=1)
#        zfaces = np.zeros_like(xmid)
#        for i in range(0, zfaces.shape[0]):
#            zfaces[i] = self.__call__([xmid[i], ymid[i]])
#
#        plt.tripcolor(x, y, triangles, facecolors=zfaces) #, edgecolors='k')



#    def plot(self, nlevel=10, vmin=None, vmax=None):
#        if vmin is None:
#            vmin = self.control.min()
#        if vmax is None:
#            vmax = self.control.max()
#        levels = np.linspace(vmin, vmax, nlevel)
#        for T_id in range(0, self.n_patchs):
#            points   = self.points[T_id,...]
#
#    #        for i in range(0, points.shape[0]):
#    #            plt.plot(points[i,0], points[i,1], "or")
#            triangles = self.triangles[T_id, ...]
##            plt.triplot(points[:,0], points[:,1], triangles, 'b-')
#
#            z = []
#            for i in range(0, points.shape[0]):
#                z.append(self.__call__(points[i,:]))
#            z = np.array(z)
#            x = points[:,0]
#            y = points[:,1]
#            plt.tricontourf(x, y, triangles, z)#, levels=levels)
#
##            z = []
##            for i in range(0, triangles.shape[0]):
##                A = points[i,0,:]
##                B = points[i,1,:]
##                C = points[i,2,:]
##
##                M = (A + B + C )/3.
##                v = self.__call_(M)
##                z.append(v)
##            z = np.array(z)
##            plt.tripcolor(x, y, triangles, facecolors=z) #, edgecolors='k')
