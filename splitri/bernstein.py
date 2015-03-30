# -*- coding: UTF-8 -*-
import numpy as np
from scipy.misc import factorial

class bernstein(object):
    def __init__(self, degree):
        self._degree = degree

        _arr = np.array(range(0, self.degree+1))
        self._factorial = factorial(_arr, exact=False)

    @property
    def degree(self):
        return self._degree

    def __call__(self, ij, x_barycentric):
        """
        evaluates the ijk bernstein polynomial at the point x
        given by its barycentric coordinates
        """
        x_barycentric = np.array(x_barycentric)
        k = self.degree - ij[0] - ij[1]
#        print "ijk ", ij, k
        w = 1. - x_barycentric.sum()

        B = x_barycentric[0]**ij[0] \
          * x_barycentric[1]**ij[1] \
          * w**k

        B *= self._factorial[self.degree]
        B /= self._factorial[ij[0]]
        B /= self._factorial[ij[1]]
        B /= self._factorial[k]

        return B
