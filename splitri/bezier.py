# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
from matplotlib.tri import TriRefiner, UniformTriRefiner
import matplotlib.pyplot as plt

from .bernstein import *
from .utils.triangle import barycentric_coords
from .trirefine import UniformBezierTriRefiner
from .utils.utils import iround


class Bezier(object):
    """
    An abstract Bezier object class.

    This Bezier class allows for the definition of Bezier surfaces or NURBS
    on triangles. First, the logical (B-net) is constructed. Then, the user can
    specify the control points and associated values for the Bernstein
    representation. If triangles are not specified, a Delaunay triangulation is
    constructed.

    Args:
        degree (int): degree of the Bernstein representation

        x (array_like):  x coordinates for vertices

        y (array_like):  y coordinates for vertices

        triangles (array_like, optional): list of 3-tuples defining the vertices indices


    Create a quarter circle NURBS curve with 2D control points and
    rational weigths and check error:

    >>> C = [[0, 1], [1, 1], [1, 0]] # 3x2 grid of 2D control points
    >>> w = [1, np.sqrt(2)/2, 1]     # rational weigths
    >>> U = [0,0,0, 1,1,1]           # knot vector

    """
    def __init__(self, degree, x, y, triangles=None):
        """
        Create a Bezier object.
        """
        self._x         = x
        self._y         = y
        self._triangles = triangles
        self._degree    = degree
        self._control_x = None
        self._control_y = None

        self._create_bezier_patch()

        # initialize control points to domain points
        self._control_x   = self.triang_ref.x.copy()
        self._control_y   = self.triang_ref.y.copy()

        self._b_coeff = np.zeros_like(self.triang_ref.x)

        self._bernstein = bernstein(self.degree)

    @property
    def b_coeff(self):
        """
        return the Bernstein coefficients for points domain
        """
        return self._b_coeff

    def control(self, i_vertex):
        """
        return the control point as a numpy array, for a given domain point
        """
        return np.array([self._control_x[i_vertex], self._control_y[i_vertex]])

    def set_control(self, i_vertex, x, y):
        """
        Set the control point by giving its coordinates x,y for a given domain point
        """
        self._control_x[i_vertex] = x
        self._control_y[i_vertex] = y

    def _create_bezier_patch(self):
        triang = tri.Triangulation(self.x, self.y, self.triangles)
        refiner = UniformBezierTriRefiner(triang)
        triang_ref, ancestors = refiner.refine_triangulation(degree=self.degree,
                                                         ancestors=True)

        self._triang = triang
        self._triang_ref = triang_ref
        self._ancestors = ancestors
        self._refiner = refiner

    def find_vertex_domain(self, P):
        """
        returns the vertex id of a given domain point.
        research is performed on the B-net (refined triangulations)
        """
        x_ref = self.triang_ref.x
        y_ref = self.triang_ref.y

        list_i = [i for i in range(0, x_ref.shape[0]) \
                  if x_ref[i]==P[0] and y_ref[i]==P[1]]
        return list_i

    def find_vertex(self, P):
        """
        returns the vertex id of a vertex.
        research is performed on the initial triangulation
        """
        x_ref = self.triang.x
        y_ref = self.triang.y

        list_i = [i for i in range(0, x_ref.shape[0]) \
                  if x_ref[i]==P[0] and y_ref[i]==P[1]]
        return list_i

    def find_simplex(self, P, tol=1.e-7):
        """
        returns the id of the triangle that contains the point P.
        If not found, returns None
        """
        ntri = self.triang.triangles.shape[0]
        vertices = np.zeros((3,2))
        for i in range(0, ntri):
            vertices[:,0] = self.triang.x[self.triang.triangles[i]]
            vertices[:,1] = self.triang.y[self.triang.triangles[i]]
            c = barycentric_coords(vertices, P)
            if (c[0]+tol >= 0.) and (c[1]+tol >= 0.) and (c[0]+c[1] <= 1.+tol):
                return i
        return None

#    def _compute_ij_ref(self, T_ref_id):
#        """
#        computes the local index of a domain-point within the initial triangulation
#        """
#        d = self.degree
#
#        triangle_ref = self.triang_ref.triangles[T_ref_id]
#        ancestor     = self.ancestors[T_ref_id]
#        triangle     = self.triang.triangles[ancestor]
#
#        x     = self.triang.x[triangle]         ; y     = self.triang.y[triangle]
#        x_ref = self.triang_ref.x[triangle_ref] ; y_ref = self.triang_ref.y[triangle_ref]
#        a11 = x[0] - x[2]
#        a21 = x[1] - x[2]
#        a12 = y[0] - y[2]
#        a22 = y[1] - y[2]
#        delta = a11*a22-a12*a21
#        X = x_ref - x[2]  ; Y = y_ref - y[2]
#        i = - X*a22 + Y*a12 ; j = - X*a21 + Y*a11
#        i /= delta        ; j /= delta
#        i *= d            ; j *= d
#        i = iround(i)     ; j = iround(j)
#
#        return np.array(i),np.array(j)

    def _compute_ij(self, T_id):
        """
        computes the local index of a domain-point within the initial triangulation
        """
        d = self.degree

        mask_ref = np.where(self.ancestors == T_id)
        triangles_ref = self.triang_ref.triangles[mask_ref]

        x   = self.triang.x[self.triang.triangles[T_id]]
        y   = self.triang.y[self.triang.triangles[T_id]]

        a11 = x[0] - x[2] ; a21 = y[0] - y[2]
        a12 = x[1] - x[2] ; a22 = y[1] - y[2]

        delta = a11*a22-a12*a21

        list_i = [] ; list_j = []
        for T_ref_id in range(0, len(triangles_ref)):
            x_ref = self.triang_ref.x[triangles_ref[T_ref_id]]
            y_ref = self.triang_ref.y[triangles_ref[T_ref_id]]

            X = x_ref - x[2]  ; Y = y_ref - y[2]

            i = X*a22 - Y*a12 ; j = -( X*a21 - Y*a11)
            i /= delta        ; j /= delta
            i *= d            ; j *= d
            ii = i            ; jj = j
            i = iround(i)     ; j = iround(j)

            ll_condition = np.all(i+j<=self.degree)
            ll_condition = ll_condition and np.all(i>=0)
            ll_condition = ll_condition and np.all(j>=0)
#            if not ll_condition:
#                print self.degree
#                print ii,jj, i, j, i+j
#                print "x_ref, y_ref ", x_ref, y_ref
#                print "x,y ", x, y
#                plt.plot(x,y,"ob")
#                plt.plot(x_ref,y_ref,"or")
#                plt.show()

            assert(ll_condition)

            list_i.append(i) ; list_j.append(j)

        return np.array(list_i), np.array(list_j)

    def evaluate_on_triangle(self, x_bary, T_id):
        """
        evaluates the Bezier surface on the triangle T_id at the point x_bary
        given by its barycentric coordinates with respect to triangle T_id
        """
        triangle = self.triang.triangles[T_id]
        mask_ref = np.where(self.ancestors == T_id)
#        print mask_ref
        triangles_ref = self.triang_ref.triangles[mask_ref]

        list_x_ref = self.triang_ref.x[triangles_ref]
        list_y_ref = self.triang_ref.y[triangles_ref]

        list_cx_ref = self._control_x[triangles_ref]
        list_cy_ref = self._control_y[triangles_ref]

        list_coeff = self.b_coeff[triangles_ref]

        list_i, list_j = self._compute_ij(T_id)
#        print list_i, list_j

        x_size = 3*len(list_x_ref)
        A = np.zeros((x_size,2))
        A[:,0] = np.array(list_x_ref.reshape(x_size))
        A[:,1] = np.array(list_y_ref.reshape(x_size))
        a = np.ascontiguousarray(A)
        unique_a, idx = np.unique(a.view([('', a.dtype)]*a.shape[1]), return_index=True)

        control = np.zeros((x_size,2))
        control[:,0] = np.array(list_cx_ref.reshape(x_size))
        control[:,1] = np.array(list_cy_ref.reshape(x_size))

        IJ = np.zeros((x_size,2), dtype=np.int32)
        IJ[:,0] = np.array(list_i.reshape(x_size))
        IJ[:,1] = np.array(list_j.reshape(x_size))

        x_ref  = np.array(list_x_ref.reshape(x_size))
        y_ref  = np.array(list_y_ref.reshape(x_size))
        coeffs = np.array(list_coeff.reshape(x_size))

        IJ = IJ[idx,:]
        coeffs = coeffs[idx]
        control = control[idx]

        x_ref  = x_ref[idx]
        y_ref  = y_ref[idx]
        vertices = np.zeros((x_ref.shape[0],2))
#        print x_ref
#        print y_ref
        vertices[:,0] = x_ref
        vertices[:,1] = y_ref

#        sum_IJ = IJ[:,0] + IJ[:,1]
#        print "============"
#        print IJ[:,0]
#        print IJ[:,1]
#        print sum_IJ[:]
#        assert np.all(sum_IJ >= 0)

        value    = 0.
        position = np.zeros(2)
        for ij in range(0, IJ.shape[0]):
            i = IJ[ij,0] ; j = IJ[ij,1]

            bern = self._bernstein([i,j], x_bary)

            value    += bern * coeffs[ij]
            position += bern * control[ij,:]

        # test if the position is the current triangle
#        if (self.find_simplex(position) == T_id):
#            print ">>> warning: coordinates ", x_bary, " out of triangle ", T_id
#            print "<<< computed position is ", position

        return value, position

    @property
    def triang(self):
        """
        returns the initial triangulation
        """
        return self._triang

    @property
    def refiner(self):
        """
        returns the triangulation refiner object
        """
        return self._refiner

    @property
    def triang_ref(self):
        """
        returns the refined triangulation (B-net)
        """
        return self._triang_ref

    @property
    def ancestors(self):
        """
        returns ancestor for each small triangle of the B-net within the initial
        triangulation
        """
        return self._ancestors

    @property
    def triangles(self):
        """
        returns initial triangles if given
        """
        return self._triangles

    @property
    def degree(self):
        """
        returns the total polynomial degree
        """
        return self._degree

    @property
    def x(self):
        """
        returns the x coordinates of vertices
        """
        return self._x

    @property
    def y(self):
        """
        returns the y coordinates of vertices
        """
        return self._y
