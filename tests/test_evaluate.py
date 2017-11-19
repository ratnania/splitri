# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from splitri.gallery.examples import square,collela \
        , domain_1, domain_2,annulus, two_triangles \
        , hexa_meshes, hexa_meshes_annulus, hexa_meshes_2
from splitri.bezier import Bezier
from splitri.utils.triangle import barycentric_coords
from matplotlib.tri import Triangulation, UniformTriRefiner
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from igakit.cad import circle
from caid.cad_geometry import cad_geometry, cad_nurbs

def test_1():
    L = 2.
    n = 10
    degree = 3

#    triang = collela(n)
    triang = square(n)

    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    x = Bzr.triang_ref.x
    y = Bzr.triang_ref.y

    Bzr._b_coeff = x**2+y**2

#    npts = 100
#    pts = np.random.random((npts,2))
#    pts = 2*pts-1.

#    old_settings = np.seterr(all='ignore')  #seterr to known value
#    np.seterr(over='raise')

    from splitri.boxsplines import triangulation_I
    n_I = 13 ; degree_I = 3 # degree
    tri_I = triangulation_I(n_I,degree_I,xmin=-1.,xmax=1.)

    triang_I = tri_I.triang
#    triang_I = tri_I.triang_ref
    npts = triang_I.x.shape[0]

    positions_x = []
    positions_y = []
    values    = []
    for i in range(0, npts):
        P = np.array([triang_I.x[i], triang_I.y[i]])

        # find in which triangle it belongs
        T_id = Bzr.find_simplex(P)
        if T_id is not None:
            # compute barycentric coordinates
            vertices = np.zeros((3,2))
            vertices[:,0] = Bzr.triang.x[Bzr.triang.triangles[T_id]]
            vertices[:,1] = Bzr.triang.y[Bzr.triang.triangles[T_id]]

            C = barycentric_coords(vertices, P)

            # evaluate the bezier surface on the point C within the triangle T
            v,position = Bzr.evaluate_on_triangle(C, T_id)
            positions_x.append(position[0])
            positions_y.append(position[1])
            values.append(v)

#    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
#    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")

    # create a triangulation for the computed positions
#    triang_new = tri.Triangulation(positions_x,positions_y)
    triang_new = triang_I
    plt.triplot(triang_new, '-', lw=0.5, color="green")

    refiner = UniformTriRefiner(triang_new)
    triang_new, values= refiner.refine_field(values, subdiv=3)

    values = np.array(values)
    print(values.min(), values.max())

    plt.tripcolor(triang_new, values, shading='gouraud', cmap=plt.cm.rainbow)


#    xlim = [-1.2,1.2] ; ylim = [-1.2,1.2]
#    plt.xlim(*xlim)  ; plt.ylim(*ylim)
    plt.show()

def test_2():
    degree = 3

    triang = two_triangles()

    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    x = Bzr.triang_ref.x
    y = Bzr.triang_ref.y

    Bzr._b_coeff = x**2+y**2

    refiner = UniformTriRefiner(Bzr.triang)
    triang_new = refiner.refine_triangulation(subdiv=3)
    npts = triang_new.x.shape[0]

    positions_x = []
    positions_y = []
    values    = []
    for i in range(0, npts):
        P = np.array([triang_new.x[i], triang_new.y[i]])

        # find in which triangle it belongs
        T_id = Bzr.find_simplex(P)
        if T_id is not None:
            # compute barycentric coordinates
            vertices = np.zeros((3,2))
            vertices[:,0] = Bzr.triang.x[Bzr.triang.triangles[T_id]]
            vertices[:,1] = Bzr.triang.y[Bzr.triang.triangles[T_id]]

            C = barycentric_coords(vertices, P)

            # evaluate the bezier surface on the point C within the triangle T
            v,position = Bzr.evaluate_on_triangle(C, T_id)
            positions_x.append(position[0])
            positions_y.append(position[1])
            values.append(v)

#    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
#    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")

    # create a triangulation for the computed positions
    triang_new = tri.Triangulation(positions_x,positions_y)
    plt.triplot(triang_new, '-', lw=0.5, color="green")

    values = np.array(values)
    print(values.min(), values.max())
    plt.tripcolor(triang_new, values, shading='gouraud', cmap=plt.cm.rainbow)

#    xlim = [-1.2,1.2] ; ylim = [-1.2,1.2]
#    plt.xlim(*xlim)  ; plt.ylim(*ylim)
    plt.show()

def test_3():
    degree = 3
    n_angles = 6
    radius = 1.
    center = np.array([0.,0.])
    angle = 2*pi/n_angles

#    triang = hexa_meshes(radius=radius, center=center, n_angles=n_angles)
    triang = hexa_meshes_annulus(min_radius=1., max_radius=2., center=center, n_angles=n_angles)

    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")

#    mask = Bzr.triang_ref.neighbors == -1
#    triangles_bnd = Bzr.triang_ref.triangles[mask]
##    print triangles_bnd
#    ancestors = Bzr.ancestors[triangles_bnd]
##    print ancestors


    nrb = circle(radius=radius, center=center, angle=angle)
    nrb_new = nrb.clone().elevate(0, times=degree-2)
    list_crv = [nrb_new]
    for i in range(1,n_angles):
        nrb = circle(radius=radius, center=center, angle=angle)
        nrb_new = nrb.clone().elevate(0, times=degree-2)
        nrb_new.rotate(i*angle)
        list_crv.append(nrb_new)

#    for nrb, col in zip(list_crv, ["b","k","y","g","r","m"]):
#        plt.plot(nrb.points[:,0],nrb.points[:,1],"-o"+col)

    for T_id in range(0, Bzr.triang.triangles.shape[0]):
        triangle = Bzr.triang.triangles[T_id]
        x = np.array(Bzr.triang.x[triangle])
        y = np.array(Bzr.triang.y[triangle])

        nrb = list_crv[T_id]

        # domain-points on the boundary are given by i=0
        # don't take into account the starting and ending points
        i = 0
        for j in range(1, Bzr.degree):
            k = Bzr.degree - j
            ijk = np.array([i,j,k], dtype=np.int32)
#            print ijk
            Px = np.dot(ijk,x) / Bzr.degree
            Py = np.dot(ijk,y) / Bzr.degree

            P = np.array([Px,Py])
            i_vertex = Bzr.find_domain_point(P)
#            print i_vertex
            Bzr.triang_ref.x[i_vertex] = nrb.points[j,0]
            Bzr.triang_ref.y[i_vertex] = nrb.points[j,1]


    x = Bzr.triang_ref.x
    y = Bzr.triang_ref.y

    Bzr._b_coeff = x**2+y**2

    refiner = UniformTriRefiner(Bzr.triang)
    triang_new = refiner.refine_triangulation(subdiv=3)
    npts = triang_new.x.shape[0]

    positions_x = []
    positions_y = []
    values    = []
    for i in range(0, npts):
        P = np.array([triang_new.x[i], triang_new.y[i]])

        # find in which triangle it belongs
        T_id = Bzr.find_simplex(P)
        if T_id is not None:
            # compute barycentric coordinates
            vertices = np.zeros((3,2))
            vertices[:,0] = Bzr.triang.x[Bzr.triang.triangles[T_id]]
            vertices[:,1] = Bzr.triang.y[Bzr.triang.triangles[T_id]]

            C = barycentric_coords(vertices, P)

            # evaluate the bezier surface on the point C within the triangle T
            v,position = Bzr.evaluate_on_triangle(C, T_id)
            positions_x.append(position[0])
            positions_y.append(position[1])
            values.append(v)

    # create a triangulation for the computed positions
    triang_new = tri.Triangulation(positions_x,positions_y)
#    plt.triplot(triang_new, '-', lw=0.5, color="green")

    values = np.array(values)
    print(values.min(), values.max())
    plt.tripcolor(triang_new, values, shading='gouraud', cmap=plt.cm.rainbow)

    xlim = [-1.2,1.2] ; ylim = [-1.2,1.2]
    plt.xlim(*xlim)  ; plt.ylim(*ylim)
    plt.show()


def test_4():
    degree = 3
    min_radius = 1.
    max_radius = 2.
    center = np.array([0.,0.])

    triang = hexa_meshes_annulus(min_radius=min_radius, max_radius=max_radius, center=center)
    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")

#    xlim = [-1.2,1.2] ; ylim = [-1.2,1.2]
#    plt.xlim(*xlim)  ; plt.ylim(*ylim)
    plt.show()


def test_5():
    degree = 3
    min_radius = 1.
    max_radius = 1.
    center = np.array([0.,0.])
    n_angles = 6
    angle = 2*pi/n_angles
    n_levels = 5

    triang = hexa_meshes_2(radius=max_radius, center=center, n_levels=n_levels)
    Bzr = Bezier(degree, triang.x, triang.y, triang.triangles)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")

    xlim = [-1.2,1.2] ; ylim = [-1.2,1.2]
    plt.xlim(*xlim)  ; plt.ylim(*ylim)
#    plt.show()

    def create_boundary(i, n_levels):
        nrb = circle(radius=max_radius, center=center, angle=angle)
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
        for enum,t in enumerate(list_t[:-1]):
            geo.split(enum,t,axis)
        geo.split(-1,list_t[-1],axis, normalize=[True,True])

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

    for T_id in range(0, len(triangles)):
        triangle = triangles[T_id]
        x = np.array(triang.x[triangle])
        y = np.array(triang.y[triangle])

        xmid = x.mean()
        ymid = y.mean()

        # associate triangle to boundary
        nrb = nearest_boundary(xmid, ymid)

#        plt.plot(nrb.points[:,0],nrb.points[:,1],"ob")
#        plt.plot(x,y,"or")


        # domain-points on the boundary are given by i=0
        # don't take into account the starting and ending points
        j = 0
        for i in range(0, Bzr.degree+1):
            k = Bzr.degree - j - i
            ijk = np.array([i,j,k], dtype=np.int32)
#            print ijk
#            print x
            Px = np.dot(ijk,x) / Bzr.degree
            Py = np.dot(ijk,y) / Bzr.degree
#            print Px.shape

            P = np.array([Px,Py])
            i_vertex = Bzr.find_domain_point(P)
#            print i_vertex
            if len(i_vertex) >0:
                i_vertex = i_vertex[0]
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
    print(values.min(), values.max())
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
