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

# ..................................
def test_1():

    radius=1.
    center=None
    n_levels=4
    degree = 3

    Bzr =  hexagonal(degree, radius=radius, center=center, n_levels=n_levels)

    for level,col in zip(range(0,n_levels+1), ["r","b", "g", "y","m","k"]):
        x,y, list_i_vertices = Bzr.get_level_vertices(level, indices=True)
#        print "level ",level," vertices found ", list_i_vertices
        plt.plot(x,y,'o'+col)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.axis("equal")
    plt.savefig("test_hexa_1.png")
# ..................................

# ..................................
def get_support(Bzr, radius, mult, center_loc=None):
    if center_loc is None:
        center_loc = Bzr.center
    Bzr_loc =  hexagonal(1, radius=mult*radius, center=center_loc, n_levels=mult)

    return Bzr_loc
# ..................................

# ..................................
def plot_support(Bzr, color=None, value=1.0):
    values = value*np.ones(len(Bzr.x))
    plt.triplot(Bzr.triang, '-', lw=0.75, color=color)
#    plt.tripcolor(Bzr.triang, values, shading='gouraud', cmap=plt.cm.rainbow)
# ..................................

# ..................................
def plot_triangle(Bzr, T_id, color=None, value=1.0):

    x_T = Bzr.triang.x[Bzr.triangles[T_id]]
    y_T = Bzr.triang.y[Bzr.triangles[T_id]]

    x_mid = 1./3. * (x_T.sum())
    y_mid = 1./3. * (y_T.sum())

    plt.plot(x_T,y_T,"or")
    plt.plot(x_mid,y_mid,"ob")
# ..................................

# ..................................
def get_non_vanishing_boxsplines(Bzr, T_id, mult, tol=1.e-7):

    # ...
    x_T = Bzr.triang.x[Bzr.triangles[T_id]]
    y_T = Bzr.triang.y[Bzr.triangles[T_id]]
    x_mid = 1./3. * (x_T.sum())
    y_mid = 1./3. * (y_T.sum())
    P = np.array([x_mid, y_mid])
    print P
    # ...

    # ... create a list of supports
    radius = Bzr.radius / Bzr.n_levels
    list_Bzr = []
    for i in range(0, len(Bzr.x)):
        center = np.array([Bzr.x[i], Bzr.y[i]])
        Bzr_loc = get_support(Bzr, radius, mult, center_loc=center)

        list_Bzr.append(Bzr_loc)
    # ...
    print len(list_Bzr)

    # ...
    list_non_vanishing_Bzr = []
    list_indices = []
    ind = 0
    for Bzr_loc in list_Bzr:
        ret = Bzr_loc.find_simplex(P, tol=tol)
        if ret is not None:
            list_non_vanishing_Bzr.append(Bzr_loc)
            ind += 1
        list_indices.append(ind)
    # ...
    print len(list_non_vanishing_Bzr)

    return list_Bzr, list_non_vanishing_Bzr, list_indices
# ..................................

# ..................................
def my_plot(Bzr, T_id, list_Bzr, ax=None):
    if ax is None:
        pt = plt
    else:
        pt = ax

    x_T = Bzr.triang.x[Bzr.triangles[T_id]]
    y_T = Bzr.triang.y[Bzr.triangles[T_id]]

    x_mid = 1./3. * (x_T.sum())
    y_mid = 1./3. * (y_T.sum())

    pt.plot(x_T,y_T,"or")
    pt.plot(x_mid,y_mid,"ob")

    pt.triplot(Bzr.triang, '-', lw=0.75, color="gray")

    for bzr in list_Bzr:
        values = np.ones(len(bzr.x))
        pt.triplot(bzr.triang, '-', lw=0.75, color="blue")
        pt.tripcolor(bzr.triang, values, shading='gouraud', cmap=plt.cm.rainbow)

# ..................................

# ..................................
def test_2():

    radius=1.
    center=None
    n_levels = 6
    mult = 3
    degree = 3 * mult - 2

    # ...
    Bzr =  hexagonal(degree, radius=radius, center=center, n_levels=n_levels)

#    plt.triplot(Bzr.triang, '-', lw=0.75, color="gray")
#    plt.triplot(Bzr.triang_ref, lw=0.5, color="green")
    # ...

    # ...
    T_id = 0

    P = np.array([0.,0.])
    P += 1.e-1
    T_id = Bzr.find_simplex(P, tol=1.e-7)
#    plot_triangle(Bzr, T_id)

    list_Bzr, list_non_vanishing_Bzr, list_indices = \
            get_non_vanishing_boxsplines(Bzr, T_id, mult, tol=1.e-7)

    print ">>> nen : ", len(list_non_vanishing_Bzr)

#    for Bzr_loc in list_non_vanishing_Bzr[:]:
#        plot_support(Bzr_loc, color="blue")
    # ...

#    plt.axis("equal")
#    plt.savefig("test_hexa_2.png")

    plt.clf()
    for enum,bzr in enumerate(list_Bzr):
        plt.clf()
        ax = plt.subplot(111)

        my_plot(Bzr, T_id, [bzr], ax=ax)

        ax.set_xlim([-1.5,1.5])
        ax.set_ylim([-1.5,1.5])
        plt.title("$nen = "+str(list_indices[enum])+"$")

#        plt.axis("equal")
        plt.axis("off")
        plt.savefig("support_"+str(enum)+".png")
# ..................................

#######################################################
if __name__ == "__main__":
#    test_1()
    test_2()
#    test_3()
#    test_4()
#    test_5()
#    test_6()
