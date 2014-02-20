import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')

fig = plt.figure(1, (7,7))
#plt.subplots_adjust(hspace=0.15,wspace=0.0)  
nx=101
ny=101

xleft=0.10
ylower = 0.08
width=0.3
height=0.3
hspace = 0.32
wspace = 0.32

colorbar_height = 0.27
colorbar_width = 0.01
# the position of the lower colorbar
colorbar_x = xleft + wspace + width + 0.01 
colorbar_y = ylower + (height-colorbar_height)/2 




# figure 1 data
fname = "p_0.dat"
p = np.loadtxt(fname,usecols=[2],unpack=True)
p=p*60
pmin = min(p)
pmax = max(p)
p = p.reshape((nx,ny))
p = np.transpose(p)


# figure 2 data
fname2 ="p_210.dat"
p2 = np.loadtxt(fname2,usecols=[2],unpack=True)
p2 = p2.reshape((nx,ny))
p2 = np.transpose(p2)
pmin0 = np.min(p2)
pmax0 = np.max(p2)
p2 = p2.reshape((nx,ny))
p2 = np.transpose(p2)

#figure 1 plot



ax1=fig.add_axes([xleft, ylower+hspace, width, height])
plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')
im1=plt.imshow(p, origin='lower',cmap='PuRd', extent=[1,nx,1,ny])
#plt.xlabel('x',labelpad=5,fontsize=20,fontproperties=tnri)
#plt.ylabel('y',labelpad=5,fontsize=20,fontproperties=tnri)
plt.xlim(1,nx)
plt.ylim(1,ny)

plt.xticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
plt.yticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
plt.text(0.15,0.85,'(a)', color='black',fontsize=16, transform=ax1.transAxes, \
        fontproperties=tnr, ha='center',va='center')
plt.text(0.8,0.8,r'$\times$60', color='black',fontsize=16, transform=ax1.transAxes, \
        fontproperties=tnr, ha='center',va='center')        
        
# upper colorbar
#ax = fig.add_axes([colorbar_x, colorbar_y+hspace, colorbar_width, colorbar_height])
#plt.tick_params(direction='out', bottom='off', top='on', left='off', right='off')
#ticks_min = np.around(pmin,decimals=5)
#ticks_max = np.around(pmax,decimals=5)
#ticks_range = np.linspace(ticks_min,ticks_max,5)
#print ticks_range
#plt.clim(ticks_min,ticks_max)
#cbar = plt.colorbar( cax=ax, ticks=ticks_range)


#cl=cbar.ax.set_yticklabels([r'2.0$\times$$10^{-5}$',r'7.8$\times$$10^{-5}$',\
                     #r'1.4$\times$$10^{-4}$',r'1.9$\times$$10^{-4}$',\
                     #r'2.5$\times$$10^{-4}$'],fontsize=13,fontproperties=tnr)


# colorbar
ax = fig.add_axes([colorbar_x, ylower, colorbar_width, height+hspace])
plt.tick_params(direction='out', bottom='off', top='on', left='off', right='off')
ticks_min = np.around(pmin0,decimals=5)
ticks_max = np.around(pmax0,decimals=5)
ticks_range = np.linspace(ticks_min,ticks_max,5)
print ticks_range
plt.clim(ticks_min,ticks_max)
cbar = plt.colorbar( cax=ax, ticks=ticks_range)

cl=cbar.ax.set_yticklabels([r'0',r'3.7$\times$$10^{-3}$',\
                     r'7.4$\times$$10^{-3}$',r'1.1$\times$$10^{-2}$',\
                     r'1.5$\times$$10^{-2}$'],fontsize=13,fontproperties=tnr)


#figure 2



ax2=fig.add_axes([xleft+wspace, ylower+hspace, width, height])
plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')
im2=plt.imshow(p2, origin='lower',cmap='PuRd', extent=[1,nx,1,ny])
plt.clim(pmin0,pmax0)
#plt.xlabel('x',labelpad=5,fontsize=20,fontproperties=tnri)
#plt.ylabel('y',labelpad=5,fontsize=20,fontproperties=tnri)
plt.xlim(1,nx)
plt.ylim(1,ny)
plt.xticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
plt.yticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
plt.text(0.15,0.85,'(b)', color='black',fontsize=16, transform=ax2.transAxes, \
        fontproperties=tnr, ha='center',va='center')

# lower colorbar
#ax = fig.add_axes([colorbar_x, colorbar_y, colorbar_width, colorbar_height])
#plt.tick_params(direction='out', bottom='off', top='on', left='off', right='off')
#ticks_min = np.around(pmin0,decimals=5)
#ticks_max = np.around(pmax0,decimals=5)
#ticks_range = np.linspace(ticks_min,ticks_max,5)
#print ticks_range
#plt.clim(ticks_min,ticks_max)
#cbar = plt.colorbar( cax=ax, ticks=ticks_range)

#cl=cbar.ax.set_yticklabels([r'0',r'3.7$\times$$10^{-3}$',\
                     #r'7.4$\times$$10^{-3}$',r'1.1$\times$$10^{-2}$',\
                     #r'1.5$\times$$10^{-2}$'],fontsize=13,fontproperties=tnr)




#figure 3
fname3 ="p_159.dat"
p3 = np.loadtxt(fname3,usecols=[2],unpack=True)
p3 = p3.reshape((nx,ny))
p3 = np.transpose(p3)


ax3=fig.add_axes([xleft, ylower, width, height])
plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')
im3=plt.imshow(p3, origin='lower',cmap='PuRd', extent=[1,nx,1,ny])
plt.clim(pmin0,pmax0)
plt.xlabel('x',labelpad=7,fontsize=20,fontproperties=tnri)
plt.ylabel('y',labelpad=2,fontsize=20,fontproperties=tnri)
plt.xlim(1,nx)
plt.ylim(1,ny)
plt.xticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
plt.yticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
plt.text(0.15,0.85,'(c)', color='black',fontsize=16, transform=ax3.transAxes, \
        fontproperties=tnr, ha='center',va='center')





#figure 4
ax4=fig.add_axes([xleft+wspace, ylower, width, height])
fname4 ="p_126.dat"
p4 = np.loadtxt(fname4,usecols=[2],unpack=True)
p4 = p4.reshape((nx,ny))
p4 = np.transpose(p4)

plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')
im4=plt.imshow(p4, origin='lower',cmap='PuRd', extent=[1,nx,1,ny])
plt.clim(pmin0,pmax0)
#plt.xlabel('x',labelpad=5,fontsize=20,fontproperties=tnri)
#plt.ylabel('y',labelpad=5,fontsize=20,fontproperties=tnri)
plt.xlim(1,nx)
plt.ylim(1,ny)
plt.xticks([1,25,50,75,101],['','','','',''], fontsize=14,fontproperties=tnr)
plt.yticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
plt.text(0.15,0.85,'(d)', color='black',fontsize=16, transform=ax4.transAxes, \
        fontproperties=tnr, ha='center',va='center')





#plt.savefig('figure3.png', bbox_inches='tight', dpi=300)
plt.show()
