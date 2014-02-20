import scipy
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import mpl
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

import matplotlib.font_manager as fm

# for defining new color scheme
import matplotlib.colors as col
import matplotlib.cm as cm


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')


width_left=0.3
width_right = 0.35
height=0.36

xmin = 0.1
ymin = 0.1

left_hspace = 0.15
LR_wspace = 0.15
subplot_hspace = 0.05
subplot_wspace = 0.0
AC_label_x_position = 0.1
AC_label_y_position = 0.7
#left_hspace_2 = 0.3


#figure A
gsA=gridspec.GridSpec(1,1)
left = xmin
right = xmin + width_left
bottom = ymin + height + left_hspace
top = bottom + height 
gsA.update(left=left,right=right, bottom=bottom, top=top )
axA = plt.subplot(gsA[0,0])
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
imA=plt.imshow(p, cmap='PuRd', origin='lower', extent=[1,nx,1,ny],aspect='auto') 
plt.clim(ticks_min,ticks_max*0.1)
#ticks_range = np.linspace(ticks_min,ticks_max*0.1,5)
print "ticks_min = " , ticks_min
print "0.1*ticks_max = " , 0.1*ticks_max
print "ticks_max = " , ticks_max
plt.xlim(1,nx)
plt.ylim(1,ny)
plt.xticks(np.linspace(1,nx,5),['0','1.25','2.50','3.75','5.00'], fontsize=16,fontproperties=tnr)
plt.yticks(np.linspace(1,ny,5), ['1', '250', '500', '750', '1001'],fontsize=16,fontproperties=tnr)
axA.set_xlabel(r'Time ($m$s)',fontsize=18,fontproperties=tnr)
axA.set_ylabel('Site index',fontsize=18,fontproperties=tnr)
axA.text(AC_label_x_position,AC_label_y_position,'(a)',fontsize=18,color='black',transform=axA.transAxes,fontproperties=tnr)




# figure C
gsC=gridspec.GridSpec(1,1)
left = xmin
right = xmin + width_left
bottom = ymin 
top = bottom + height
gsC.update(left=left,right=right, bottom=bottom, top=top )
axC = plt.subplot(gsC[0,0])
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
print "ticks_min = " , ticks_min
print "0.1*ticks_max = " , 0.3*ticks_max
print "ticks_max = " , ticks_max
plt.imshow(p, cmap='PuRd', origin='lower', extent=[1,nx,1,ninterval],aspect='auto') 
plt.clim(ticks_min,ticks_max*0.3)
ticks_range = np.linspace(ticks_min,ticks_max*0.3,5)
plt.xlim(1,nx)
plt.ylim(1,ninterval)
plt.xticks(np.linspace(1,nx,5),['0','0.45','0.90','1.35','1.80'], fontsize=16,fontproperties=tnr)
#plt.yticks(np.linspace(1,ny,5), ['300', '400', '500', '600', '700'],fontsize=12,fontproperties=tnr)
plt.yticks(np.linspace(1,ninterval,5), ['200', '350', '500', '650', '800'],fontsize=16,fontproperties=tnr)
axC.set_xlabel(r'Time ($m$s)',fontsize=18,fontproperties=tnr)
axC.set_ylabel('Site index',fontsize=18,fontproperties=tnr)
axC.text(AC_label_x_position,AC_label_y_position,'(c)',fontsize=18,color='black',transform=axC.transAxes,fontproperties=tnr)




# figure B 
gsB=gridspec.GridSpec(2,2, height_ratios=[1,5], width_ratios=[15,1])
left = xmin + width_left + LR_wspace 
right = left + width_right #+ 0.15
bottom = ymin + height + left_hspace
top = bottom + height
gsB.update(hspace=0.0, wspace=0.1, left=left,right=right, bottom=bottom, top=top )

axB = plt.subplot(gsB[0:2,0])
axB.tick_params(direction='out',bottom='on', top='off', left='on', right='off')

site0, probability0 = np.loadtxt('planewave_focusing_t_zero.txt').T
site1, probability1 = np.loadtxt('planewave_focusing_t_most_narrow_focus_501.txt').T
site2, probability2 = np.loadtxt('planewave_focusing_t_most_narrow_focus_350.txt').T

axB.set_xlim(1, 1001)
axB.set_ylim(0, 0.144)
ytick_position = np.linspace(0,0.144,5)

xtick_number=6
xtick_position =np.linspace(0,1001,xtick_number)
plt.xticks(xtick_position, ['1', '200', '400','600','800','1001'],fontsize=16,fontproperties=tnr)
plt.yticks(ytick_position, ['0', '0.036','0.072','0.108','0.144'], fontsize=16,fontproperties=tnr)           
axB.plot(site0,probability0*20,'k--', site1, probability1, 'r-', site2, probability2, 'b-')
axB.set_xlabel('Site index', fontsize=18,fontproperties=tnr)
axB.set_ylabel('Probability', fontsize=18,fontproperties=tnr)
axB.text(0.18,0.85,'(b)',fontsize=18,transform=axB.transAxes,fontproperties=tnr)
axB.text(0.09,0.16,r'$\times$ 20',fontsize=16,transform=axB.transAxes,fontproperties=tnr)

# two colorbars for figure A
axA_colorbar1 = plt.subplot(gsB[1,1])
axA_colorbar1.tick_params(direction='out',bottom='off', top='off', left='off', right='on')
cb1 = mpl.colorbar.ColorbarBase(axA_colorbar1, ticks=[0,0.5,1.0], cmap=cm.get_cmap('PuRd'))#, orientation='horizontal')
cb1.ax.set_yticklabels(['0.0','0.007',\
                     '0.014'],fontsize=13,fontproperties=tnr)
                     

axA_colorbar2 = plt.subplot(gsB[0,1])
axA_colorbar2.tick_params(direction='out',bottom='off', top='off', left='off', right='on')
cmap = mpl.colors.ListedColormap([mpl.colors.rgb2hex([0.40392157435417175, 0.0, 0.12156862765550613]), \
                                  mpl.colors.rgb2hex([0.40392157435417175, 0.0, 0.12156862765550613])])
cb1 = mpl.colorbar.ColorbarBase(axA_colorbar2, cmap=cmap, ticks=[0,1.0]) #, orientation='horizontal')
cb1.ax.set_yticklabels(['','0.138'],fontsize=13,fontproperties=tnr)


#figure D
gsD=gridspec.GridSpec(2,2, height_ratios=[1,5], width_ratios=[15,1])
left = xmin + width_left + LR_wspace 
right = left + width_right #+ 0.15
bottom = ymin 
top = bottom + height 
gsD.update(hspace=0.0, wspace=0.1, left=left,right=right, bottom=bottom, top=top )
axD = plt.subplot(gsD[0:2,0])
axD.tick_params(direction='out',bottom='on', top='off', left='on', right='off')
site0, probability0 = np.loadtxt('wavepacket_focusing_t_zero.txt').T
site1, probability1 = np.loadtxt('wavepacket_focusing_t_most_narrow_focus_501.txt').T
site2, probability2 = np.loadtxt('wavepacket_focusing_t_most_narrow_focus_600.txt').T
axD.set_xlim(1, 1001)
axD.set_ylim(0, 0.036)
ytick_position = np.linspace(0,0.036,5)
xtick_number=6
xtick_position =np.linspace(0,1001,xtick_number)
plt.xticks(xtick_position, ['1', '200', '400','600','800','1001'],fontsize=16,fontproperties=tnr)
plt.yticks(ytick_position, ['0', '0.009','0.018','0.027','0.036'], fontsize=16,fontproperties=tnr)           
axD.plot(site0,probability0*5,'k--', site1, probability1, 'r-', site2, probability2, 'b-')
axD.set_xlabel('Site index', fontsize=18,fontproperties=tnr)
axD.set_ylabel('Probability', fontsize=18,fontproperties=tnr)
axD.text(0.18,0.85,'(d)',fontsize=18,transform=axD.transAxes,fontproperties=tnr)
axD.text(0.09,0.14,r'$\times$ 5',fontsize=16,transform=axD.transAxes,fontproperties=tnr)

# two colorbar for figure C
axC_colorbar1 = plt.subplot(gsD[1,1])
axC_colorbar1.tick_params(direction='out',bottom='off', top='off', left='off', right='on')
cb1 = mpl.colorbar.ColorbarBase(axC_colorbar1, ticks=[0,0.5,1.0], cmap=cm.get_cmap('PuRd') ) #), orientation='horizontal')
cb1.ax.set_yticklabels(['0.0','0.005',\
                     '0.011'],fontsize=13,fontproperties=tnr)

axC_colorbar2 = plt.subplot(gsD[0,1])
axC_colorbar2.tick_params(direction='out',bottom='off', top='off', left='off', right='on')
cmap = mpl.colors.ListedColormap([mpl.colors.rgb2hex([0.40392157435417175, 0.0, 0.12156862765550613]), \
                                  mpl.colors.rgb2hex([0.40392157435417175, 0.0, 0.12156862765550613])])
cb1 = mpl.colorbar.ColorbarBase(axC_colorbar2, cmap=cmap, ticks=[0,1.0])  #, orientation='horizontal')
cb1.ax.set_yticklabels(['','0.035'],fontsize=13,fontproperties=tnr)





plt.savefig('focusing-1d.pdf', bbox_inches='tight',dpi=300)
plt.show()

