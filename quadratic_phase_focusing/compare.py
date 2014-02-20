import matplotlib.pyplot as plt
import numpy as np

nx = 101
ny = 101
#x1, y1, p = np.loadtxt("t_250.dat", unpack=True)

plt.figure(figsize=(7,12))

         
plt.subplots_adjust(hspace=0.2)
ax=plt.subplot(211)

plt.text(0.5, 1.2, 'Maximum probability enhanced by a factor of 56', \
         color='red', \
         horizontalalignment='center', \
         verticalalignment='center', \
         transform = ax.transAxes)

ntime = 0
fname = "p_" + str(ntime) + ".dat"
p = np.loadtxt(fname,usecols=[2],unpack=True)
pmin = min(p)
pmax = max(p)
#print pmax
p = p.reshape((nx,ny))
p = np.transpose(p)
plt.title(r'Time = 0 $ms$',fontsize=18)	
plt.imshow(p, origin='lower',extent=[1,nx,1,ny])

ticks_min = np.around(pmin,decimals=5)
ticks_max = np.around(pmax,decimals=5)
plt.clim(ticks_min,ticks_max)
ticks_range = np.linspace(ticks_min,ticks_max,5)

plt.colorbar(ticks=ticks_range)
#plt.colorbar()
plt.xlabel('x',fontsize=18)
plt.ylabel('y',fontsize=18)
plt.xlim(1,nx)
plt.ylim(1,ny)
plt.xticks(np.linspace(1,nx,5))
plt.yticks(np.linspace(1,ny,5))

plt.subplot(212)
ntime = 159
fname = "p_" + str(ntime) + ".dat"
p = np.loadtxt(fname,usecols=[2],unpack=True)
pmin = min(p)
pmax = max(p)
#print pmax
p = p.reshape((nx,ny))
p = np.transpose(p)
plt.title(r'Time = 1.05 $ms$',fontsize=18)	
plt.imshow(p, origin='lower',extent=[1,nx,1,ny])

ticks_min = np.around(pmin,decimals=5)
ticks_max = np.around(pmax,decimals=5)
plt.clim(ticks_min,ticks_max)
ticks_range = np.linspace(ticks_min,ticks_max,5)

plt.colorbar(ticks=ticks_range)
#plt.colorbar()
plt.xlabel('x',fontsize=18)
plt.ylabel('y',fontsize=18)
plt.xlim(1,nx)
plt.ylim(1,ny)
plt.xticks(np.linspace(1,nx,5))
plt.yticks(np.linspace(1,ny,5))
plt.savefig('focusing_1ms_sigma_1d-3_focus_60.60.pdf')
plt.show()


