import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import os

def format_e(n):
	if abs(n-0)<1.e-10:
		return '0.00'
			
	a = '%.2e' % n 
	return str(a)

#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')



def makeMovie(directory, timePoints, saveTo):
	for t in timePoints:
		filename = directory + "p_" + str(t) + ".dat"
		#time = "time = " + str(t) + " $\mu s$"
		nx=101
		ny=101
		p = np.loadtxt(filename,usecols=[2],unpack=True)
		pmin = 0.0 #min(p)
		pmax = max(p)
		p = p.reshape((nx,ny))
		p = np.transpose(p)
		
		
		plt.figure(t)
		plt.tick_params(direction='out', bottom='off', top='off', left='off', right='off')
		plt.imshow(p, origin='lower',cmap='PuRd', extent=[1,nx,1,ny])
		#plt.xlabel('x',labelpad=10,fontsize=20,fontproperties=tnri)
		#plt.ylabel('y',labelpad=5,fontsize=20,fontproperties=tnri)
		plt.xlim(1,nx)
		plt.ylim(1,ny)
		#plt.set_xticks([])
		#plt.set_yticks([])
		plt.xticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
		plt.yticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
		
		# colorbar
		#label_min = np.around(pmin*1000,decimals=1)
		#label_max = np.around(pmax*1000,decimals=1)
		ticks_range = np.linspace(pmin,pmax,4)
		ticks_label=map(format_e, ticks_range)
		
		plt.clim(pmin,pmax)
		cbar = plt.colorbar( ticks=ticks_range)
		cbar.ax.set_yticklabels(ticks_label,fontsize=13,fontproperties=tnr)
	
		#plt.show()
		if 0<=t<=9:
			fig = saveTo +"00" + str(t) + ".png"
		elif 10<=t<=99:
			fig = saveTo +"0" + str(t) + ".png"
		else:
			fig = saveTo + str(t) + ".png"
	
		plt.savefig(fig)
		plt.close()
	
	make_movie = "convert -delay 20 " + saveTo+"*.png " + saveTo+".gif" 
	os.system(make_movie)
	os.system("rm -rf *.png")



dir1='/home/pxiang/run_jobs/wavepacket_1d/wavepacket_focusing_2d/wavepacket_focusing_Ja_22.83_kHz/time_0.5ms_sigma_1d-3_focus_51.51/'
gifFile="toCenter"
maxTime = 293
timePoints = range(0,51)+ range(50,293,10)
makeMovie(dir1,timePoints, gifFile)

dir2='/home/pxiang/run_jobs/wavepacket_1d/wavepacket_focusing_2d/wavepacket_focusing_Ja_22.83_kHz/time_0.5ms_sigma_1d-3_focus_31.31/'
gifFile="toLowerLeft"
timePoints = range(0,51)+ range(50,270,12)
makeMovie(dir2,timePoints,gifFile)
