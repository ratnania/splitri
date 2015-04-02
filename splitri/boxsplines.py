# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import scipy
from trirefine import UniformBezierTriRefiner
from splitri.triangulation import triangulation_square_I, triangulation_square_II

class BoxSpline(object):
    def __init__(self, triang):
        """
        triang is a triangulation I or II
        """

        self._b_coeff = np.zeros_like(triang.triang_ref.x)

    @property
    def b_coeff(self):
        return self._b_coeff


class BoxSpline_211(BoxSpline):
    def __init__(self, triang, center=[0.,0.]):
        if not isinstance(triang, triangulation_square_I):
            raise ValueError("Expected a Triangulation of type I object")

        assert (triang.degree==2)
        BoxSpline.__init__(self, triang)

        dx = triang.dx / triang.degree  ; dy = triang.dy / triang.degree

        center = np.asarray(center) + np.array([dx,0.])
        list_P = [center]                         ;  list_v = [2.]
        list_P.append(center+np.array([dx,0.]))   ;  list_v.append(1.)
        list_P.append(center+np.array([0.,dy]))   ;  list_v.append(1.)
        list_P.append(center+np.array([dx,dy]))   ;  list_v.append(1.)
        list_P.append(center+np.array([-dx,0.]))  ;  list_v.append(1.)
        list_P.append(center+np.array([0.,-dy]))  ;  list_v.append(1.)
        list_P.append(center+np.array([-dx,-dy])) ;  list_v.append(1.)

        for P,v in zip(list_P, list_v):
            i = triang.find_domain_point(P)
            self._b_coeff[i] = v
        self._b_coeff /= 2.

class BoxSpline_221(BoxSpline):
    def __init__(self, triang, center=[0.,0.]):
        if not isinstance(triang, triangulation_square_I):
            raise ValueError("Expected a Triangulation of type I object")

        assert (triang.degree==3)
        BoxSpline.__init__(self, triang)

        dx = triang.dx / triang.degree  ; dy = triang.dy / triang.degree

        center = np.asarray(center)

        list_P = []                               ;  list_v = []
        # 1st column
        i = -1
        for j,v in zip( [0, 1, 2] \
                       ,[1.,1.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 2nd column
        i = 0
        for j,v in zip( [-1, 0, 1, 2, 3] \
                       ,[1.,2.,3.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 3rd column
        i = 1
        for j,v in zip( [-1, 0, 1, 2, 3, 4] \
                       ,[1.,3.,4.,4.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 4th column
        i = 2
        for j,v in zip( [-1, 0, 1, 2, 3, 4] \
                       ,[1.,2.,4.,4.,3.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 5th column
        i = 3
        for j,v in zip( [0, 1, 2, 3, 4] \
                       ,[1.,2.,3.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 6th column
        i = 4
        for j,v in zip( [1, 2, 3] \
                       ,[1.,1.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)


        for P,v in zip(list_P, list_v):
            i = triang.find_domain_point(P)
            self._b_coeff[i] = v
        self._b_coeff /= 6.


class BoxSpline_222(BoxSpline):
    def __init__(self, triang, center=[0.,0.]):
        if not isinstance(triang, triangulation_square_I):
            raise ValueError("Expected a Triangulation of type I object")

        assert (triang.degree==4)
        BoxSpline.__init__(self, triang)

        dx = triang.dx / triang.degree  ; dy = triang.dy / triang.degree

        center = np.asarray(center)

        list_P = []                               ;  list_v = []
        #
        i = -1
        for j,v in zip( [0, 1, 2, 3] \
                       ,[1.,1.,1.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 0
        for j,v in zip( [-1, 0, 1, 2, 3, 4, 5] \
                       ,[1.,2.,3.,4.,3.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 1
        for j,v in zip( [-1, 0, 1, 2, 3, 4, 5, 6] \
                       ,[1.,3.,4.,6.,6.,4.,3.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 2
        for j,v in zip( [-1, 0, 1, 2,  3, 4, 5, 6, 7] \
                       ,[1.,4.,6.,8.,10.,8.,6.,4.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 3
        for j,v in zip( [-1, 0, 1,  2,  3,  4,  5, 6, 7, 8] \
                       ,[1.,3.,6.,10.,12.,12.,10.,6.,4.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 4
        for j,v in zip( [0, 1,  2,  3,  4,  5, 6, 7, 8] \
                       ,[2.,4.,8.,12.,12.,12.,8.,4.,2.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 5
        for j,v in zip( [0, 1,  2,  3,  4,  5, 6, 7, 8, 9] \
                       ,[1.,3.,6.,10.,12.,12.,10.,6.,3.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 6
        for j,v in zip( [1, 2,  3, 4, 5, 6, 7, 8, 9] \
                       ,[1.,4.,6.,8.,10.,8.,6.,4.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 7
        for j,v in zip( [2, 3, 4, 5, 6, 7, 8, 9] \
                       ,[1.,3.,4.,6.,6.,4.,3.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 8
        for j,v in zip( [3, 4, 5, 6, 7, 8, 9] \
                       ,[1.,2.,3.,4.,3.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        #
        i = 9
        for j,v in zip( [5, 6, 7, 8] \
                       ,[1.,1.,1.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)

        for P,v in zip(list_P, list_v):
            i = triang.find_domain_point(P)
            self._b_coeff[i] = v
        self._b_coeff /= 24.


class BoxSpline_1111(BoxSpline):
    def __init__(self, triang, center=[0.,0.]):
        if not isinstance(triang, triangulation_square_II):
            raise ValueError("Expected a Triangulation of type II object")

        assert (triang.degree==2)
        BoxSpline.__init__(self, triang)

        N = (2*triang.degree)

        dx = triang.dx / N ; dy = triang.dy / N

        center = np.asarray(center)

        list_P = []                               ;  list_v = []
        # 1st column
        i = -1
        for j,v in zip( [ 1, 2] \
                       ,[1.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 2nd column
        i = 0
        for j,v in zip( [0, 1, 2, 3] \
                       ,[1.,2.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 3rd column
        i = 1
        for j,v in zip( [0, 1, 2, 3] \
                       ,[1.,2.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 4th column
        i = 2
        for j,v in zip( [ 1, 2] \
                       ,[1.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)

        for P,v in zip(list_P, list_v):
            i = triang.find_domain_point(P)
            self._b_coeff[i] = v
        # TODO scaling to be corrected
        self._b_coeff /= 2.


class BoxSpline_2111(BoxSpline):
    def __init__(self, triang, center=[0.,0.]):
        if not isinstance(triang, triangulation_square_II):
            raise ValueError("Expected a Triangulation of type II object")

        assert (triang.degree==3)
        BoxSpline.__init__(self, triang)

        N = 2*triang.degree + 1

        dx = triang.dx / N ; dy = triang.dy / N

        center = np.asarray(center)

        list_P = []                               ;  list_v = []
        # 1st column
        i = -1
        for j,v in zip( [ 1, 2] \
                       ,[1.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 2nd column
        i = 0
        for j,v in zip( [0, 1, 2, 3] \
                       ,[1.,2.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 3rd column
        i = 1
        for j,v in zip( [0, 1, 2, 3] \
                       ,[1.,2.,2.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)
        # 4th column
        i = 2
        for j,v in zip( [ 1, 2] \
                       ,[1.,1.]):
            list_P.append(center+np.array([i*dx,j*dy])) ;  list_v.append(v)

        for P,v in zip(list_P, list_v):
            i = triang.find_domain_point(P)
            self._b_coeff[i] = v

        self._b_coeff /= 48.


