# coding: utf-8
#import splitri.core.boxsplines as bx
#import numpy as np
#
#my_size=501
#bsize1=0; bsize2=0
#b = np.zeros((my_size,my_size), dtype=np.int32)
#a = np.zeros((my_size,my_size), dtype=np.int32)
#
##l=2 ; m=1; n=1
##a,b,bsize1,bsize2=bx.bnet.bs3dm(l,m,n,my_size)
#
#k=2; l=2 ; m=1; n=1
#a,b,bsize1,bsize2=bx.bnet.bs4dm(k,l,m,n,my_size)
#b = b[0:bsize1+1, 0:bsize2+1]
#import matplotlib.pyplot as plt
#plt.contourf(b)
#plt.colorbar()
#plt.show()


from splitri.boxsplines import Bnet

bnet1 = Bnet(3,3,3)
print bnet1.coeff
bnet2 = Bnet(2,2,2,2)
print bnet2.coeff

import matplotlib.pyplot as plt
plt.pcolor(bnet1.coeff)
plt.colorbar()
plt.show()
