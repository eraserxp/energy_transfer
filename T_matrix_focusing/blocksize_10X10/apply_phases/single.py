import matplotlib.pyplot as plt
import numpy as np

nx = 101
ny = 101
#x1, y1, p = np.loadtxt("t_250.dat", unpack=True)
ntime = 1

while ntime:
  ntime = raw_input("ntime = ")
  fname = "p_" + ntime + ".dat"
  p = np.loadtxt(fname,usecols=[2],unpack=True)
  pmin = min(p)
  pmax = max(p)
  print pmax
  p = p.reshape((nx,ny))
  p = np.transpose(p)
#  plt.title('t=4.e-7 s')	
  plt.imshow(p, origin='lower',extent=[1,nx,1,ny])

  ticks_min = np.around(pmin,decimals=5)
  ticks_max = np.around(pmax,decimals=5)
  plt.clim(ticks_min,ticks_max)
  ticks_range = np.linspace(ticks_min,ticks_max,5)

  plt.colorbar(ticks=ticks_range)
  #plt.colorbar()
  plt.xlabel('x')
  plt.ylabel('y')
  plt.xlim(1,nx)
  plt.ylim(1,ny)
  plt.xticks(np.linspace(1,nx,5))
  plt.yticks(np.linspace(1,ny,5))

  plt.show()
#plt.savefig('plot1.pdf')

