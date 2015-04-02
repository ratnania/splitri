# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import scipy
from .trirefine import UniformBezierTriRefiner
from .bezier import Bezier

__author__="ARA"
__all__ = ['triangulation_square_I', 'triangulation_square_II', 'hexagonal']

class square_boxsplines(object):
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

    def find_domain_point(self, P):
        x_ref = self.triang_ref.x
        y_ref = self.triang_ref.y

        list_i = [i for i in range(0, x_ref.shape[0]) \
                  if x_ref[i]==P[0] and y_ref[i]==P[1]]
        return list_i


class triangulation_square_I(square_boxsplines):
    def __init__(self,n,degree,xmin=-1.,xmax=1.):
        square_boxsplines.__init__(self,n,degree,xmin=xmin,xmax=xmax)

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

class triangulation_square_II(square_boxsplines):
    def __init__(self,n,degree,xmin=-1.,xmax=1.):
        square_boxsplines.__init__(self,n,degree,xmin=xmin,xmax=xmax)

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


class hexagonal(Bezier):
    """
    create an hexagonal mesh object
    """
    def __init__(self, degree, radius=1., center=None, n_levels=2):

        n_angles = 6
        if center is None:
            center = np.array([0.,0.])

        self._radius = radius
        self._center = center
        self._n_levels = n_levels

        angles = np.linspace(0.,2*np.pi,n_angles+1)
        # construct points
        x = radius * cos(angles)
        y = radius * sin(angles)

        points = np.ones((n_angles+1,2))
        points *= center
        points[1:,0] += x[:-1]
        points[1:,1] += y[:-1]
        # construct triangles
        I0 = 0 # center
        triangles = []
        for i in range(0, n_angles):
            I1 = i+1
            I2 = I1+1
            if I2==n_angles+1:
                I2 = 1
            triangles.append([I0,I1,I2])

        Bzr = Bezier(n_levels, points[:,0], points[:,1], triangles)
        triang = Bzr.triang_ref

        Bezier.__init__(self, degree, triang.x, triang.y, triang.triangles)

    @property
    def radius(self):
        return self._radius

    @property
    def center(self):
        return self._center

    @property
    def n_levels(self):
        return self._n_levels

    def get_level_vertices(self, level, indices=False):
        """
        returns vertices on the level number level.
        the returned array has a size of
        6*(level+1)
        returns indices of these points on the triangulation, if indices is True
        """
        assert(level <= self.n_levels)

        if level == 0:
            if not indices:
                return [self.center[0]],  [self.center[1]]
            else:
                i_vertex  = self.find_vertex(self.center)
                return [self.center[0]],  [self.center[1]],i_vertex

        radius = level * self.radius / (self.n_levels)
        n_angles = 6
        angles = np.linspace(0.,2*np.pi,n_angles+1)

        # construct points
        x = radius * cos(angles)
        y = radius * sin(angles)

        x_new = [] ; y_new = [] ; list_i_vertices = []
        for i in range(0, n_angles):
            # remove the last point
            t = np.linspace(0.,1.,level+1)[0:-1]

            x_tmp = (1.-t) * x[i] + t * x[i+1]
            y_tmp = (1.-t) * y[i] + t * y[i+1]

            i_vertices = np.zeros(x_tmp.shape[0],dtype=np.int32)
            for j in range(0, x_tmp.shape[0]):
                P = np.array([x_tmp[j],y_tmp[j]])
                i_vertex  = self.find_vertex(P)
                plt.plot(P[0],P[1],"dk")
                if len(i_vertex) > 0:
                    i_vertices[j] = i_vertex[0]

            x_new.append(x_tmp)
            y_new.append(y_tmp)
            list_i_vertices.append(i_vertices)

        x_new = np.array(x_new)
        y_new = np.array(y_new)
        i_vertices= np.array(list_i_vertices, dtype=np.int32)

        if not indices:
            return x_new,y_new
        else:
            return x_new,y_new,i_vertices

    def get_label_domain_points(self, level, label, indices=False):
        """
        returns domain points on the label level.
        the returned array has a size of
        6*(level+1)
        returns indices of these points on the triangulation, if indices is True
        """

        if level == 0:
            if not indices:
                return [self.center[0]],  [self.center[1]]
            else:
                i_vertex  = self.find_vertex(self.center)
                return [self.center[0]],  [self.center[1]],i_vertex

        n_label_levels = self.n_levels * self.degree
        i_label_level  = (level-1) * self.degree + label

        assert(i_label_level <= n_label_levels)

        radius = i_label_level * self.radius / n_label_levels
        n_angles = 6
        angles = np.linspace(0.,2*np.pi,n_angles+1)

        # construct points
        x = radius * cos(angles)
        y = radius * sin(angles)

        x_new = [] ; y_new = [] ; list_i_vertices = []
        for i in range(0, n_angles):
            # remove the last point
            t = np.linspace(0.,1.,i_label_level+1)[0:-1]

            x_tmp = (1.-t) * x[i] + t * x[i+1]
            y_tmp = (1.-t) * y[i] + t * y[i+1]

            i_vertices = np.zeros(x_tmp.shape[0],dtype=np.int32)
            for j in range(0, x_tmp.shape[0]):
                P = np.array([x_tmp[j],y_tmp[j]])
                i_vertex  = self.find_domain_point(P)
                plt.plot(P[0],P[1],"dk")
                if len(i_vertex) > 0:
                    i_vertices[j] = i_vertex[0]

            x_new.append(x_tmp)
            y_new.append(y_tmp)
            list_i_vertices.append(i_vertices)

        x_new = np.array(x_new)
        y_new = np.array(y_new)
        i_vertices= np.array(list_i_vertices, dtype=np.int32)

        if not indices:
            return x_new,y_new
        else:
            return x_new,y_new,i_vertices
