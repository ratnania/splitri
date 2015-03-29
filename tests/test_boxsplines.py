# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from splitri.boxsplines import *

def test_211():
    L = 2.
    n = 9

    e = [2,1,1] ; degree = 2

    plot_field = True

    tri_box = triangulation_I(n, degree)
    Box     = BoxSpline_211(tri_box)

    if plot_field:
        str_e = ""
        for _e in e:
            str_e+=str(_e)
        title = "Box-Splines $B_{"+str_e+"}$"

    #    plt.tripcolor(triang_ref, B_coeff, shading='gouraud', cmap=plt.cm.rainbow)

        refiner = UniformTriRefiner(tri_box.triang_ref)
        tri_refi, z_test_refi = refiner.refine_field(Box.b_coeff, subdiv=4)
        plt.tripcolor(tri_refi, z_test_refi, shading='gouraud', cmap=plt.cm.rainbow) ; plt.colorbar()

        plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
        plt.triplot(tri_box.triang_ref, lw=0.5, color="white")
        plt.title(title)
        plt.show()
    else:
        plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
        plt.triplot(tri_box.triang_ref, lw=0.5, color="blue")
        plt.show()


def test_221():
    L = 2.
    n = 9

    e = [2,2,1] ; degree = 3

    plot_field = True

    tri_box = triangulation_I(n, degree)
    Box     = BoxSpline_221(tri_box)

    if plot_field:
        str_e = ""
        for _e in e:
            str_e+=str(_e)
        title = "Box-Splines $B_{"+str_e+"}$"

    #    plt.tripcolor(triang_ref, B_coeff, shading='gouraud', cmap=plt.cm.rainbow)

        refiner = UniformTriRefiner(tri_box.triang_ref)
        tri_refi, z_test_refi = refiner.refine_field(Box.b_coeff, subdiv=4)
        plt.tripcolor(tri_refi, z_test_refi, shading='gouraud', cmap=plt.cm.rainbow) ; plt.colorbar()

        plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
        plt.triplot(tri_box.triang_ref, lw=0.5, color="white")
        plt.title(title)
        plt.show()
    else:
        plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
        plt.triplot(tri_box.triang_ref, lw=0.5, color="blue")
        plt.show()

def test_222():
    L = 2.
    n = 9

    e = [2,2,2] ; degree = 4

    plot_field = True

    tri_box = triangulation_I(n, degree)
    Box     = BoxSpline_222(tri_box)

    if plot_field:
        str_e = ""
        for _e in e:
            str_e+=str(_e)
        title = "Box-Splines $B_{"+str_e+"}$"

    #    plt.tripcolor(triang_ref, B_coeff, shading='gouraud', cmap=plt.cm.rainbow)

        refiner = UniformTriRefiner(tri_box.triang_ref)
        tri_refi, z_test_refi = refiner.refine_field(Box.b_coeff, subdiv=4)
        plt.tripcolor(tri_refi, z_test_refi, shading='gouraud', cmap=plt.cm.rainbow) ; plt.colorbar()

        plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
        plt.triplot(tri_box.triang_ref, lw=0.5, color="white")
        plt.title(title)
        plt.show()
    else:
        plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
        plt.triplot(tri_box.triang_ref, lw=0.5, color="blue")
        plt.show()

def test_1111():
    L = 2.
    n = 9

    e = [1,1,1,1] ; degree = 2

    plot_field = True

    tri_box = triangulation_II(n, degree)
    Box     = BoxSpline_1111(tri_box)

    if plot_field:
        str_e = ""
        for _e in e:
            str_e+=str(_e)
        title = "Box-Splines $B_{"+str_e+"}$"

    #    plt.tripcolor(triang_ref, B_coeff, shading='gouraud', cmap=plt.cm.rainbow)

        refiner = UniformTriRefiner(tri_box.triang_ref)
        tri_refi, z_test_refi = refiner.refine_field(Box.b_coeff, subdiv=4)
        plt.tripcolor(tri_refi, z_test_refi, shading='gouraud', cmap=plt.cm.rainbow) ; plt.colorbar()

        plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
        plt.triplot(tri_box.triang_ref, lw=0.5, color="white")
        plt.title(title)
        plt.show()
    else:
        plt.triplot(tri_box.triang, '-', lw=0.75, color="red")
        plt.triplot(tri_box.triang_ref, lw=0.5, color="blue")
        plt.show()

#######################################################
if __name__ == "__main__":
    test_211()
    test_221()
    test_222()
    test_1111()
