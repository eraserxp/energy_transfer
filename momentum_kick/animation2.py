import matplotlib.pyplot as plt
import numpy as np
import os

for t in np.arange(0,1001,10):
	filename = "x_" + str(t) + ".dat"
	filename2 = "k_" + str(t) + ".dat"
	time = "time = " + str(t*0.003) + " $\mu s$"

	x,p = np.loadtxt(filename,unpack=True)
	k,p2 = np.loadtxt(filename2,unpack=True)
#	pmin = min(p)
#	pmax = max(p)
	plt.figure(t)
	plt.subplots_adjust(hspace=0.3)
	
	plt.subplot(211)
	plt.plot(k,p2,'blue')
	plt.title(time)
	plt.xlabel('ka')
	plt.ylabel('Probability')
	plt.xlim(-np.pi,np.pi)
	plt.ylim(0,0.06)
	plt.text(1, 0.05, 'k space', bbox=dict(facecolor='blue', alpha=0.5))
		

	plt.subplot(212)
	plt.plot(x,p,'r')
	plt.xlabel('Molecule index')
	plt.ylabel('Probability')
	plt.xlim(1,1001)
	plt.ylim(0,0.05)
	plt.text(650, 0.04, 'x space', bbox=dict(facecolor='red', alpha=0.5))	


	if 0<=t<=9:
		fig = "plot000" + str(t) + ".png"
	elif 10<=t<=99:
		fig = "plot00" + str(t) + ".png"
	elif 100<=t<=999:
		fig = "plot0"+ str(t) + ".png"
	else:
		fig = "plot" + str(t) + ".png"

	plt.savefig(fig)
	plt.close(t)

make_movie = "convert -delay 20 plot*.png momentum_kick.gif"
os.system(make_movie)
os.system("rm plot*.png")


