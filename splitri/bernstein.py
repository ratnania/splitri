# -*- coding: UTF-8 -*-
import numpy as np
from scipy.misc import factorial

class bernstein(object):
    def __init__(self, degree):
        self._degree = degree

        _arr = np.array(list(range(0, self.degree+1)))
        self._factorial = factorial(_arr, exact=False)

    @property
    def degree(self):
        return self._degree

    def __call__(self, ij, x_barycentric, debug=False):
        """
        evaluates the ijk bernstein polynomial at the point x
        given by its barycentric coordinates
        """
        x_barycentric = np.array(x_barycentric)
        i = ij[0]
        j = ij[1]
#        if i==self.degree:
#            j = 0
        k = self.degree - i - j
        assert(i>=0)
        assert(j>=0)
        assert(k>=0)

        w = 1. - x_barycentric.sum()

#        print x_barycentric,w,ij,k

        B = x_barycentric[0]**i \
          * x_barycentric[1]**j \
          * w**k

        B *= self._factorial[self.degree]
        B /= self._factorial[i]
        B /= self._factorial[j]
        B /= self._factorial[k]

        if debug:
            if not np.isfinite(B):
                print(("Problem occurs with ",  x_barycentric,w,i,j,k, B))

        return B
