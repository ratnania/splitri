# -*- coding: UTF-8 -*-
import numpy as np
from scipy.special import binom as sc_binom

class bernstein(object):
    def __init__(self, degree):
        self._degree = degree
        self._initialize_binom()

    @property
    def degree(self):
        return self._degree

    def _initialize_binom(self):
        allBinom = []
        for d in range(0,self.degree+1):
            values = np.zeros(d+1)
            for i in range(0, d+1):
                values[i] = sc_binom(d,i)
            allBinom.append(values)
        self._allBinom = allBinom

#    for d in range(0, degree+1):
#        print(len(allBinom[d]))

    def _binom(self, n,i):
        return self._allBinom[n][i]

    def __call__(self, ijk, x_barycentric):
        """
        evaluates the ijk bernstein polynomial at the point x
        given by its barycentric coordinates
        """
        v = x_barycentric[0]**ijk[0] \
          * x_barycentric[1]**ijk[1] \
          * x_barycentric[2]**ijk[2]

        v *= self._binom(self.degree,ijk[0])
        v *= self._binom(self.degree-ijk[0],ijk[1])
        return v
