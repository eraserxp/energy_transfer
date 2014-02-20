from pylab import *

#Z = rand(6,6)

#subplot(2,1,1)
#c = pcolor(Z)
#title('default: no edges')

#subplot(2,1,2)

Z=loadtxt("phase_matrix.txt", usecols=[2]).T
Z=Z.reshape(20,20)

c = pcolor(Z, edgecolors='k', linewidths=1,cmap='hsv')
#xticks([0,1,2,3,4,5,6],['1','20','40','60','80','90','101'])
#yticks([0,1,2,3,4,5,6],['1','20','40','60','80','90','101'])
colorbar()
clim(0,2*pi)
title('thick edges')

show()
