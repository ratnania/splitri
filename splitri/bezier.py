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
    def __init__(self, degree, x, y, triangles=None):

        self._x         = x
        self._y         = y
        self._triangles = triangles
        self._degree    = degree
        self._control   = None

        self._create_bezier_patch()

        self._b_coeff = np.zeros_like(self.triang_ref.x)

        for i in range(0, 5):
            self._compute_ij(i)

    @property
    def b_coeff(self):
        return self._b_coeff

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
        returns the vertex id
        research is performed on the refined
        triangulations with the domain points
        """
        x_ref = self.triang_ref.x
        y_ref = self.triang_ref.y

        list_i = [i for i in range(0, x_ref.shape[0]) \
                  if x_ref[i]==P[0] and y_ref[i]==P[1]]
        return list_i

    def find_vertex(self, P):
        """
        returns the vertex id
        research is performed on the initial triangulation
        """
        x_ref = self.triang.x
        y_ref = self.triang.y

        list_i = [i for i in range(0, x_ref.shape[0]) \
                  if x_ref[i]==P[0] and y_ref[i]==P[1]]
        return list_i

    def find_simplex(self, P, tol=1.e-7):
        """
        returns the id of the triangle that contains the point xyz
        """
        ntri = self.triang.triangles.shape[0]
        vertices = np.zeros((3,2))
        for i in range(0, ntri):
            vertices[:,0] = self.triang.x[self.triang.triangles[i]]
            vertices[:,1] = self.triang.y[self.triang.triangles[i]]
            c = barycentric_coords(vertices, P)
#            k = np.where(c>=0)
#            i = k[0].min()
            if (c[0]+tol >= 0.) and (c[1]+tol >= 0.) and (c[0]+c[1] <= 1.+tol):
                return i
        return None

    def _compute_ij(self, T_ref_id):
        """
        computes the local index of a domain-point within the initial triangulation
        """
#        print self.ancestors
#        print self.triang_ref.triangles
        d = self.degree
        triangle_ref = self.triang_ref.triangles[T_ref_id]
        ancestor     = self.ancestors[T_ref_id]
        triangle     = self.triang.triangles[ancestor]

        x     = self.triang.x[triangle]         ; y     = self.triang.y[triangle]
        x_ref = self.triang_ref.x[triangle_ref] ; y_ref = self.triang_ref.y[triangle_ref]
        a11 = x[0] - x[2]
        a21 = x[1] - x[2]
        a12 = y[0] - y[2]
        a22 = y[1] - y[2]
        delta = a11*a22-a12*a21
#        print "-----------"
#        print T_ref_id
#        print triangle_ref
#        print triangle
#        print "x", x
#        print "y", y
#        print "x_ref", x_ref
#        print "y_ref", y_ref
        X = x_ref - x[2]  ; Y = y_ref - y[2]
        i = - X*a22 + Y*a12 ; j = - X*a21 + Y*a11
        i /= delta        ; j /= delta
        i *= d            ; j *= d
        i = iround(i)     ; j = iround(j)
        print i,j
        plt.plot(x,y,"oy")
        plt.plot(x_ref,y_ref,"xr")
        return i,j
#        print "-----------"


#+    def __call__(self, xy):
#+        """
#+        evaluates the surface using bernstein polynomial
#+        """
#+        A = np.array(xy)
#+        T_id = self.find_simplex(A)
#+        vertices = self.vertices[T_id,...]
#+        coeffs = self.control[T_id,...]
#+        c = barycentric_coords(vertices, A)
#+
#+        value = 0.
#+        i_pos = 0
#+        for i in range(0, self.degree+1):
#+            for j in range(0, self.degree+1-i):
#+                k = self.degree - i - j
#+                bern = self._bernstein([i,j,k], c)
#+                value += bern * coeffs[i_pos]
#+                i_pos += 1
#+        return value

#    def evaluate_on_triangle(self, x_bary, T_id):
#        """
#        evaluates the Bezier surface on the triangle T_id at the point x_bary
#        """
#        triangle = self.triang.triangles[T_id]
#        mask_ref = np.where(self.triang_ref.triangles_ancestors == T_id)
#        triangles_ref = self.triang_ref.triangles[mask_ref]
#        value = 0.
#        i_pos = 0
#        for T_ref_id in range(0, triangles_ref.shape[0]):
#            i,j = self.get_refined_ij(T_ref_id)
#            k = self.degree - i - j
#            bern = self._bernstein([i,j,k], x_bary)
#            value += bern * coeffs[i_pos]
#            i_pos += 1
#        return value

    @property
    def triang(self):
        return self._triang

    @property
    def refiner(self):
        return self._refiner

    @property
    def triang_ref(self):
        return self._triang_ref

    @property
    def ancestors(self):
        return self._ancestors

    @property
    def triangles(self):
        return self._triangles

    @property
    def degree(self):
        return self._degree

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def control(self):
        return self.triang_ref.x, self.triang_ref.y


