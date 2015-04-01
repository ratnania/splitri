# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from splitri.utils.triangle import barycentric_coords
from splitri.triangulation import triangulation_square_I, triangulation_square_II
import splitri.triangulation as stri
from splitri.triangulation import hexagonal
from matplotlib.tri import Triangulation, UniformTriRefiner
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
    n_levels=5
    degree = 3

    Bzr =  hexagonal(degree, radius=radius, center=center, n_levels=n_levels)

    plt.triplot(Bzr.triang, '-', lw=0.75, color="red")
    plt.triplot(Bzr.triang_ref, lw=0.5, color="blue")
    plt.show()

#######################################################
if __name__ == "__main__":
    test_1()
    test_2()
    test_3()
