# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from splitri.utils.triangle import barycentric_coords
from splitri.utils.utils import construct_curve_from_points
from splitri.triangulation import triangulation_square_I, triangulation_square_II
import splitri.triangulation as stri
from splitri.triangulation import hexagonal
from matplotlib.tri import Triangulation, UniformTriRefiner
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from igakit.cad import circle
from caid.cad_geometry import cad_geometry, cad_nurbs

def test_1():
    n = 10
    degree = 3

    Bzr = triangulation_square_I(n, degree)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()

def test_2():
    n = 10
    degree = 3

    Bzr = triangulation_square_II(n, degree)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()


def test_3():

    radius=1.
    center=None
    n_levels=4
    degree = 3

    Bzr =  hexagonal(degree, radius=radius, center=center, n_levels=n_levels)

    for level,col in zip(range(0,n_levels+1), ["r","b", "g", "y","m","k"]):
        x,y, list_i_vertices = Bzr.get_level_vertices(level, indices=True)
        print "level ",level," vertices found ", list_i_vertices
        plt.plot(x,y,'o'+col)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()

def test_4():

    radius=1.
    center=None
    n_levels=4
    degree = 3

    Bzr =  hexagonal(degree, radius=radius, center=center, n_levels=n_levels)

    for level,col in zip(range(0,n_levels+1), ["r","b", "g", "y","m","k"]):
        for label in range(0, degree+1):
            x,y, list_i_vertices = Bzr.get_label_domain_points(level, label, indices=True)
            print "level ",level," vertices found ", list_i_vertices
            plt.plot(x,y,'o'+col)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()

def test_5():

    radius=1.
    center=None
    n_levels=7
    degree = 3

    Bzr =  hexagonal(degree, radius=radius, center=center, n_levels=n_levels)

#    for level,col in zip(range(0,n_levels+1), ["r","b", "g", "y","m","k"]):
#        for label in range(0, degree+1):
#            x,y, list_i_vertices = Bzr.get_label_domain_points(level, label, indices=True)
#            print "level ",level," vertices found ", list_i_vertices
#            plt.plot(x,y,'o'+col)
#
#    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
#    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
#    plt.show()

    def cloud_data(n_levels):
        """
        generates a list of list of points, each array describes 1/6 of a flux
        surface
        """
        center = np.array([2.,0.])

        n_angles = 6
        angles = np.linspace(0.,2*np.pi,n_angles+1)
        radius = 1.

        list_pts = []
        for level in range(1, n_levels+1):
            _radius = level * radius / (n_levels+1)
            sl = []
            for angle_b,angle_e in zip(angles[0:-1],angles[1:]):
                _angles = np.linspace(angle_b,angle_e,50)
                x = center[0] + _radius * np.cos(_angles)
                y = center[1] + _radius * np.sin(_angles)
                sl.append([x,y])
            list_pts.append(sl)
        return list_pts

    list_pts = cloud_data(n_levels*degree)
    pts_level = list_pts[0]
    pts_section = pts_level[0]
    pts = pts_section[0]
    print len(list_pts),len(pts_level), len(pts_section), len(pts)

#    for i in range(0,len(list_pts)):
#        pts_level = list_pts[i]
#        pts_section = pts_level[0]
#        pts = pts_section[0]
#        level = i + 1
#        for pts_section in pts_level:
#            plt.plot(pts_section[0], pts_section[1],".")
#    plt.show()

    i_label_level = 0
    for level in range(1,n_levels+1):
        for label in range(0,degree+1):
#            i_label_level = (level-1) * degree + label

            if i_label_level == n_levels*degree:
                break

            pts_level = list_pts[i_label_level]
            for pts_section in pts_level:
                x = pts_section[0]
                y = pts_section[1]
                px = 3
                knots = [0.25,0.5,0.75]
                geo = construct_curve_from_points(x,y, px, \
                                knots=knots, method='chord', alpha=0.1, nx=None)
                plt.plot(geo[0].points[:,0],geo[0].points[:,1],'-o')
                print geo[0].knots[0]
            i_label_level += 1
    plt.show()


def test_6():

    radius=1.
    center=None
    n_levels=5
    degree = 3

    center = np.array([0.,0.])
    n_angles = 6
    angle = 2*pi/n_angles
    n_levels = 5

    Bzr =  hexagonal(degree, radius=radius, center=center, n_levels=n_levels)


    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")

    xlim = [-1.2,1.2] ; ylim = [-1.2,1.2]
    plt.xlim(*xlim)  ; plt.ylim(*ylim)
#    plt.show()

    def create_boundary(i, n_levels):
        nrb = circle(radius=radius, center=center, angle=angle)
        axis = 0
        times = nrb.degree[axis]
        list_t = np.linspace(0.,1.,n_levels+2)[1:-1]
        for t in list_t:
            for _ in range(0, times):
                nrb.insert(axis,t)
        nrb = nrb.clone().elevate(axis, times=degree-2)
        if i > 0:
            nrb.rotate(i*angle)
        # split the curve into a list of bezier curves
        geo = cad_geometry()
        cad_nrb = cad_nurbs(nrb.knots, nrb.points, weights=nrb.weights)
        geo.append(cad_nrb)
#        for enum,t in enumerate(list_t[:-1]):
#            geo.split(enum,t,axis)
#        geo.split(-1,list_t[-1],axis, normalize=[True,True])

#        for nrb in geo:
#            plt.plot(nrb.points[:,0],nrb.points[:,1],"-ob")
#        plt.plot(geo[1].points[:,0],geo[1].points[:,1],"-or")
#        plt.show()
        return geo

    def nearest_boundary(xmid, ymid):
        # associate triangle to boundary
        dist = np.zeros(len(list_crv))
        for enum, nrb in enumerate(list_crv):
            dist[enum] = (xmid-nrb_centers[enum,0])**2 \
                       + (ymid-nrb_centers[enum,1])**2
        i_nrb = dist.argmin()
        return list_crv[i_nrb]

    list_crv = []
    for i in range(0,n_angles):
        geo = create_boundary(i, n_levels)
        for nrb in geo:
            list_crv.append(nrb)

#    for nrb, col in zip(list_crv, ["b","k","y","g","r","m"]):
#        plt.plot(nrb.points[:,0],nrb.points[:,1],"-o"+col)
#    plt.show()

    triang = Bzr.triang
    _mask = np.array(triang.neighbors == -1)
    mask = np.array([_m.any() for _m in _mask])
    mask_id = np.where(mask)
#    print triang.triangles
#    print mask
#    print triang.triangles[mask]
    triangles = triang.triangles[mask_id]
#    print triangles

    nrb_centers = np.zeros((len(list_crv),2))
    for i in range(0, len(list_crv)):
        nrb = list_crv[i]
        nrb_centers[i,0] = nrb.points[:,0].mean()
        nrb_centers[i,1] = nrb.points[:,1].mean()

    level = n_levels
    x,y, list_i_vertices = Bzr.get_level_vertices(level, indices=True)

    for i_level in range(0,x.shape[0]):
        _x = x[i_level]
        _y = y[i_level]
        i_vertices = list_i_vertices[i_level]

        nrb = list_crv[i_level]
        for i in range(0, _x.shape[0]):
            P = np.array([_x[i],_y[i]])
            i_vertex = i_vertices[i]
            print i, i_vertex, nrb.points.shape[0]
            Bzr.set_control(i_vertex, nrb.points[i,0], nrb.points[i,1])

#    plt.plot(Bzr.triang_ref.x, Bzr.triang_ref.y,"xr")

    x = Bzr.triang_ref.x
    y = Bzr.triang_ref.y

    Bzr._b_coeff = x**2+y**2

    refiner = UniformTriRefiner(Bzr.triang)
    subdiv=1
    triang_new = refiner.refine_triangulation(subdiv=subdiv)
    npts = triang_new.x.shape[0]

    positions_x = []
    positions_y = []
    values    = []
    for i in range(0, npts):
        P = np.array([triang_new.x[i], triang_new.y[i]])

        # find in which triangle it belongs
        T_id = Bzr.find_simplex(P)
        if T_id is not None:
#            print "Point ",P, " found in triangle ", T_id
            # compute barycentric coordinates
            vertices = np.zeros((3,2))
            vertices[:,0] = Bzr.triang.x[Bzr.triang.triangles[T_id]]
            vertices[:,1] = Bzr.triang.y[Bzr.triang.triangles[T_id]]

            C = barycentric_coords(vertices, P)
#            print "bary coord ", C
#            plt.plot(vertices[:,0],vertices[:,1],"dy")
#            plt.plot(P[0],P[1],"dg")
#            plt.show()

            # evaluate the bezier surface on the point C within the triangle T
            v,position = Bzr.evaluate_on_triangle(C, T_id)
            positions_x.append(position[0])
            positions_y.append(position[1])
            values.append(v)

    # create a triangulation for the computed positions
    triang_new = tri.Triangulation(positions_x,positions_y,triang_new.triangles)
#    plt.triplot(triang_new, '-', lw=0.5, color="green")

    values = np.array(values)
    print values.min(), values.max()
    plt.tripcolor(triang_new, values, shading='gouraud', cmap=plt.cm.rainbow)
    plt.colorbar()

    xlim = [-1.2,1.2] ; ylim = [-1.2,1.2]
    plt.xlim(*xlim)  ; plt.ylim(*ylim)
    plt.show()


#######################################################
if __name__ == "__main__":
#    test_1()
#    test_2()
#    test_3()
#    test_4()
    test_5()
#    test_6()
