# -*- coding: UTF-8 -*-
import numpy as np
import numpy.linalg as la

# ...
def barycentric_coords(vertices, point):
    """
    computes the barycentric coordinates of the point (2d) with respect to the
    triangle defined by vertices
    returns a 2d array
    the 3rd component is given by 1-v.sum()
    """
    T = (np.array(vertices[:-1])-vertices[-1]).T
    v = np.dot(la.inv(T), np.array(point)-vertices[-1])
#    v.resize(len(vertices))
#    v[-1] = 1-v.sum()
#    print v
    return v
# ...


