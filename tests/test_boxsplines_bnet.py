# coding: utf-8
import splitri.core.boxsplines as bx
import numpy as np

my_size=21
bsize1=0; bsize2=0
b = np.zeros((my_size,my_size), dtype=np.int32)
a = np.zeros((my_size,my_size), dtype=np.int32)

l=2 ; m=1; n=1
a,b,bsize1,bsize2=bx.bnet.bs3dm(l,m,n,my_size)

print b.transpose()
