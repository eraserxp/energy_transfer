import matplotlib.pyplot as plt
import numpy as np

#x1, y1, p = np.loadtxt("t_250.dat", unpack=True)
p = np.loadtxt("t_300.dat",usecols=[2],unpack=True)
pmin = min(p)
pmax = max(p)
print pmax
p = p.reshape((1001,1001))
p = np.transpose(p)
plt.title('t=4.e-7 s')	
plt.imshow(p, origin='lower',extent=[-500,500,-500,500])

ticks_min = np.around(pmin,decimals=5)
ticks_max = np.around(pmax,decimals=5)
plt.clim(ticks_min,ticks_max)
ticks_range = np.linspace(ticks_min,ticks_max,6)

#plt.colorbar(ticks=ticks_range)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-500,500)
plt.ylim(-500,500)
plt.xticks(np.arange(-500,600,100))
plt.yticks(np.arange(-500,600,100))

plt.show()
#plt.savefig('plot1.pdf')

