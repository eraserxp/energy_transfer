import scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

import matplotlib.font_manager as fm


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')


# obtain the data and produce the images: 
#'focusing_wavepacket.png' and 'focusing_planewave.png' 

#plt.figure(figsize=(3+3.0/8, 3+3.0/8))
gs=gridspec.GridSpec(2,2,width_ratios=[2,2])
#gs.update(wspace=0.37,hspace=0.36)
gs.update(wspace=0.49,hspace=0.32,left=0.16,right=0.85) #hspace=0.32



#FIGURE C
ax1=plt.subplot(gs[2])
plt.tick_params(direction='out',bottom='on', top='off', left='on', right='off')

nx = 1001
ny = 1001

fname = "wavepacket_distribution_vs_time.txt"

p2 = np.loadtxt(fname,usecols=[2],unpack=True)
p2 = p2.reshape((nx,ny))
p2 = np.transpose(p2)

# cut a part of p2
nstart = 200 #300
nend = 800 #700
ninterval = nend-nstart+1
p = p2[nstart:nend,:]
print np.shape(p)
pmin = np.min(p)
pmax = np.max(p)

ticks_min = np.around(pmin,decimals=5)
ticks_max = np.around(pmax,decimals=5)

plt.imshow(p, cmap='PuRd', origin='lower', extent=[1,nx,1,ninterval],aspect='auto') 
#plt.clim(ticks_min,ticks_max*0.1)
#ticks_range = np.linspace(ticks_min,ticks_max*0.1,5)
plt.clim(ticks_min,ticks_max*0.3)
ticks_range = np.linspace(ticks_min,ticks_max,5)
plt.colorbar(ticks=[0, ticks_max*0.1, ticks_max*0.2, ticks_max*0.3, ticks_max*0.5, ticks_max])


plt.xlim(1,nx)
#plt.ylim(1,ny)
plt.ylim(1,ninterval)
plt.xticks(np.linspace(1,nx,5),['0','0.45','0.90','1.35','1.80'], fontsize=16,fontproperties=tnr)
#plt.yticks(np.linspace(1,ny,5), ['300', '400', '500', '600', '700'],fontsize=12,fontproperties=tnr)
plt.yticks(np.linspace(1,ninterval,5), ['200', '350', '500', '650', '800'],fontsize=16,fontproperties=tnr)



###############################

ax1.set_xlabel(r'Time ($m$s)',fontsize=18,fontproperties=tnr)
ax1.set_ylabel('Site index',fontsize=18,fontproperties=tnr)
ax1.text(0.2,0.85,'(c)',fontsize=18,color='black',transform=ax1.transAxes,fontproperties=tnr)










# FIGURE D
ax2=plt.subplot(gs[3])
ax2.tick_params(direction='out',bottom='on', top='off', left='on', right='off')

site0, probability0 = np.loadtxt('wavepacket_focusing_t_zero.txt').T
site1, probability1 = np.loadtxt('wavepacket_focusing_t_most_narrow_focus_501.txt').T
site2, probability2 = np.loadtxt('wavepacket_focusing_t_most_narrow_focus_600.txt').T

ax2.set_xlim(1, 1001)
ax2.set_ylim(0, 0.036)
ytick_position = np.linspace(0,0.036,5)

xtick_number=6
xtick_position =np.linspace(0,1001,xtick_number)
plt.xticks(xtick_position, ['1', '200', '400','600','800','1001'],fontsize=16,fontproperties=tnr)
plt.yticks(ytick_position, ['0', '0.009','0.018','0.027','0.036'], fontsize=16,fontproperties=tnr)           
ax2.plot(site0,probability0*5,'k--', site1, probability1, 'r-', site2, probability2, 'b-')
ax2.set_xlabel('Site index', fontsize=18,fontproperties=tnr)
ax2.set_ylabel('Probability', fontsize=18,fontproperties=tnr)
ax2.text(0.18,0.85,'(d)',fontsize=18,transform=ax2.transAxes,fontproperties=tnr)
ax2.text(0.09,0.14,r'$\times$ 5',fontsize=16,transform=ax2.transAxes,fontproperties=tnr)





#FIGURE A
ax3=plt.subplot(gs[0])
plt.tick_params(direction='out',bottom='on', top='off', left='on', right='off')

nx = 1001
ny = 1001

fname = "planewave_distribution_vs_time.txt"

p = np.loadtxt(fname,usecols=[2],unpack=True)
p = p.reshape((nx,ny))
p = np.transpose(p)


print np.shape(p)
pmin = np.min(p)
pmax = np.max(p)

ticks_min = np.around(pmin,decimals=5)
ticks_max = np.around(pmax,decimals=5)

plt.imshow(p, cmap='PuRd', origin='lower', extent=[1,nx,1,ny],aspect='auto') 
plt.clim(ticks_min,ticks_max*0.1)
ticks_range = np.linspace(ticks_min,ticks_max*0.1,5)


plt.xlim(1,nx)
plt.ylim(1,ny)
plt.xticks(np.linspace(1,nx,5),['0','1.25','2.50','3.75','5.00'], fontsize=16,fontproperties=tnr)
plt.yticks(np.linspace(1,ny,5), ['1', '250', '500', '750', '1001'],fontsize=16,fontproperties=tnr)


###############################

ax3.set_xlabel(r'Time ($m$s)',fontsize=18,fontproperties=tnr)
ax3.set_ylabel('Site index',fontsize=18,fontproperties=tnr)
ax3.text(0.2,0.85,'(a)',fontsize=18,color='black',transform=ax3.transAxes,fontproperties=tnr)






# FIGURE B
ax4=plt.subplot(gs[1])
ax4.tick_params(direction='out',bottom='on', top='off', left='on', right='off')

site0, probability0 = np.loadtxt('planewave_focusing_t_zero.txt').T
site1, probability1 = np.loadtxt('planewave_focusing_t_most_narrow_focus_501.txt').T
site2, probability2 = np.loadtxt('planewave_focusing_t_most_narrow_focus_350.txt').T

ax4.set_xlim(1, 1001)
ax4.set_ylim(0, 0.144)
ytick_position = np.linspace(0,0.144,5)

xtick_number=6
xtick_position =np.linspace(0,1001,xtick_number)
plt.xticks(xtick_position, ['1', '200', '400','600','800','1001'],fontsize=16,fontproperties=tnr)
plt.yticks(ytick_position, ['0', '0.036','0.072','0.108','0.144'], fontsize=16,fontproperties=tnr)           
ax4.plot(site0,probability0*20,'k--', site1, probability1, 'r-', site2, probability2, 'b-')
ax4.set_xlabel('Site index', fontsize=18,fontproperties=tnr)
ax4.set_ylabel('Probability', fontsize=18,fontproperties=tnr)
ax4.text(0.18,0.85,'(b)',fontsize=18,transform=ax4.transAxes,fontproperties=tnr)
ax4.text(0.09,0.16,r'$\times$ 20',fontsize=16,transform=ax4.transAxes,fontproperties=tnr)


#plt.savefig('focusing-1d.pdf', bbox_inches='tight',dpi=300)
plt.show()
