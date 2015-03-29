# -*- coding: UTF-8 -*-
import numpy as np
from numpy import cos, sin, pi
from scipy.spatial import Delaunay
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from splitri.gallery.examples import collela
from splitri.trirefine import *

def test_1():
    L = 2.
    n = 15
    degree = 3

    triang = collela(n)

    refiner = UniformBezierTriRefiner(triang)
    triang_ref, ancestors = refiner.refine_triangulation(degree=degree,
                                                         ancestors=True)

    plt.triplot(triang, '-', lw=0.75, color="red")
    plt.triplot(triang_ref, lw=0.5, color="blue")
    plt.show()


#######################################################
if __name__ == "__main__":
    test_1()
