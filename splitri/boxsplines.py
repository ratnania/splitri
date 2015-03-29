
# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
import numpy.linalg as la
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import scipy
from trirefine import UniformBezierTriRefiner
from matplotlib.tri import TriRefiner, UniformTriRefiner

class triangulation_boxplines(object):
    def __init__(self,n,degree,xmin=-1.,xmax=1.):
        t = np.linspace(xmin,xmax,n)
        X,Y=np.meshgrid(t,t)
        x = X.reshape(X.size)
        y = Y.reshape(Y.size)

        self._dx = (xmax-xmin) / (n-1) ; self._dy = self._dx
        self._xmin = xmin
        self._xmax = xmax

        self._x = x
        self._y = y
        self._degree = degree
        self._triangles = None
        self._L = xmax-xmin

    def _create_bezier_patch(self):
        triang = tri.Triangulation(self.x, self.y, self.triangles)
        refiner = UniformBezierTriRefiner(triang)
        triang_ref, ancestors = refiner.refine_triangulation(degree=self.degree,
                                                         ancestors=True)

        self._triang = triang
        self._triang_ref = triang_ref
        self._ancestors = ancestors
        self._refiner = refiner

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
    def xmin(self):
        return self._xmin

    @property
    def xmax(self):
        return self._xmax

    @property
    def dx(self):
        return self._dx

    @property
    def dy(self):
        return self._dy

    @property
    def L(self):
        return self._L

    def find_vertex(self, P):
        x_ref = self.triang_ref.x
        y_ref = self.triang_ref.y

        list_i = [i for i in range(0, x_ref.shape[0]) \
                  if x_ref[i]==P[0] and y_ref[i]==P[1]]
        return list_i


class triangulation_I(triangulation_boxplines):
    def __init__(self,n,degree,xmin=-1.,xmax=1.):
        triangulation_boxplines.__init__(self,n,degree,xmin=xmin,xmax=xmax)

        triangles = []
        for j in range(0,n-1):
            for i in range(0,n-1):
                I1 = i+j*n ; I2 = i+1+j*n ; I3 = i+1+(j+1)*n
                T = [I1,I2,I3]
                triangles.append(T)

                I1 = i+j*n ; I2 = i+(j+1)*n ; I3 = i+1+(j+1)*n
                T = [I1,I2,I3]
                triangles.append(T)

        self._triangles = triangles

        self._create_bezier_patch()

class triangulation_II(triangulation_boxplines):
    def __init__(self,n,degree,xmin=-1.,xmax=1.):
        triangulation_boxplines.__init__(self,n,degree,xmin=xmin,xmax=xmax)

        n_mid = n-1
        t_mid = np.linspace(xmin+self.dx/2,xmax-self.dx/2, n_mid)
        X_mid,Y_mid=np.meshgrid(t_mid,t_mid)
        x_mid = X_mid.reshape(X_mid.size)
        y_mid = Y_mid.reshape(Y_mid.size)

        n_base = n * n
        x = np.concatenate((self.x,x_mid), axis=0)
        y = np.concatenate((self.y,y_mid), axis=0)

        triangles = []
        for j in range(0,n-1):
            for i in range(0,n-1):
                # triangle bottom
                I_mid = i+j*n_mid ; I1 = i+j*n ; I2 = i+1+j*n
                T = [I_mid+n_base,I1,I2]
                triangles.append(T)

                # triangle right
                I_mid = i+j*n_mid ; I1 = i+1+j*n ; I2 = i+1+(j+1)*n
                T = [I_mid+n_base,I1,I2]
                triangles.append(T)

                # triangle top
                I_mid = i+j*n_mid ; I1 = i+1+(j+1)*n ; I2 = i+(j+1)*n
                T = [I_mid+n_base,I1,I2]
                triangles.append(T)

                # triangle left
                I_mid = i+j*n_mid ; I1 = i+(j+1)*n ; I2 = i+j*n
                T = [I_mid+n_base,I1,I2]
                triangles.append(T)

        self._x = x
        self._y = y
        self._triangles = triangles

        self._create_bezier_patch()

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
    def __init__(self, triang_I, center=[0.,0.]):

        assert (triang_I.degree==2)
        BoxSpline.__init__(self, triang_I)

        dx = triang_I.dx / triang_I.degree  ; dy = triang_I.dy / triang_I.degree

        center = np.asarray(center) + np.array([dx,0.])
        list_P = [center]                         ;  list_v = [2.]
        list_P.append(center+np.array([dx,0.]))   ;  list_v.append(1.)
        list_P.append(center+np.array([0.,dy]))   ;  list_v.append(1.)
        list_P.append(center+np.array([dx,dy]))   ;  list_v.append(1.)
        list_P.append(center+np.array([-dx,0.]))  ;  list_v.append(1.)
        list_P.append(center+np.array([0.,-dy]))  ;  list_v.append(1.)
        list_P.append(center+np.array([-dx,-dy])) ;  list_v.append(1.)

        for P,v in zip(list_P, list_v):
            i = triang_I.find_vertex(P)
            self._b_coeff[i] = v
        self._b_coeff /= 2.

class BoxSpline_221(BoxSpline):
    def __init__(self, triang_I, center=[0.,0.]):

        assert (triang_I.degree==3)
        BoxSpline.__init__(self, triang_I)

        dx = triang_I.dx / triang_I.degree  ; dy = triang_I.dy / triang_I.degree

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
            i = triang_I.find_vertex(P)
            self._b_coeff[i] = v
        self._b_coeff /= 6.


class BoxSpline_222(BoxSpline):
    def __init__(self, triang_I, center=[0.,0.]):

        assert (triang_I.degree==2)
        BoxSpline.__init__(self, triang_I)

        dx = triang_I.dx / triang_I.degree  ; dy = triang_I.dy / triang_I.degree

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
            i = triang_I.find_vertex(P)
            self._b_coeff[i] = v
        self._b_coeff /= 24.
#######################################################

if __name__ == "__main__":
    L = 2.
    n = 9
    e = [2,1,1] ; degree = 2
#    e = [2,2,1] ; degree = 3
#    e = [2,2,2] ; degree = 4


    if e == [2,1,1]:
        tri_box = triangulation_I(n, degree)
        Box     = BoxSpline_211(tri_box)
    if e == [2,2,1]:
        tri_box = triangulation_I(n, degree)
        Box     = BoxSpline_221(tri_box)
    if e == [2,2,2]:
        tri_box = triangulation_I(n, degree)
        Box     = BoxSpline_222(tri_box)

    tri_box = triangulation_II(n, degree)

#    title = "Box-Splines $B_{"+str(e[0])+str(e[1])+str(e[2])+"}$"
#
##    plt.tripcolor(triang_ref, B_coeff, shading='gouraud', cmap=plt.cm.rainbow)
#
#    refiner = UniformTriRefiner(tri_box.triang_ref)
#    tri_refi, z_test_refi = refiner.refine_field(Box.b_coeff, subdiv=4)
#    plt.tripcolor(tri_refi, z_test_refi, shading='gouraud', cmap=plt.cm.rainbow) ; plt.colorbar()
#
#    plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
#    plt.triplot(tri_box.triang_ref, lw=0.5, color="white")
#    plt.title(title)
#    plt.show()

    plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
    plt.triplot(tri_box.triang_ref, lw=0.5, color="blue")
    plt.show()


