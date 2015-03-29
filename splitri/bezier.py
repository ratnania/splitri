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


class Bezier(object):
    def __init__(self, degree, x, y, triangles=None):

        self._x         = x
        self._y         = y
        self._triangles = triangles
        self._degree    = degree
        self._control   = None

        self._create_bezier_patch()

        self._b_coeff = np.zeros_like(self.triang_ref.x)

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


