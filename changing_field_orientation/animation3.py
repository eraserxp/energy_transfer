import matplotlib.pyplot as plt
import numpy as np
import os

for t in np.arange(1,301):
	filename = "t_" + str(t) + ".dat"
	time = "time = " + str(t*0.01) + " ms"

	p = np.loadtxt(filename,usecols=[2],unpack=True)
	pmin = min(p)
	pmax = max(p)
	p = p.reshape((1001,1001))
	p = np.transpose(p)
	plt.figure(t)
	plt.title(time)	
	plt.imshow(p, origin='lower',extent=[-500,500,-500,500])

	ticks_min = np.around(pmin,decimals=5)
	ticks_max = np.around(pmax,decimals=5)
	plt.clim(ticks_min,ticks_max)
	ticks_range = np.linspace(ticks_min,ticks_max,6)

	plt.colorbar(ticks=ticks_range)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.xlim(-500,500)
	plt.ylim(-500,500)
	plt.xticks(np.arange(-500,600,100))
	plt.yticks(np.arange(-500,600,100))

#	plt.show()
	if 0<t<=9:
		fig = "plot00" + str(t) + ".png"
	elif 10<=t<=99:
		fig = "plot0" + str(t) + ".png"
	else:
		fig = "plot" + str(t) + ".png"

	plt.savefig(fig)
	plt.close(t)

make_movie = "convert -delay 10 plot*.png kax_0.5pi_kay_0.5pi_theta_0_90_phi_0.gif"
os.system(make_movie)
#os.system("rm plot*.png")


