import os
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
hspace = 0.32 #0.38
wspace = 0.32 #0.42

colorbar_height = 0.27
colorbar_width = 0.01
# the position of the lower colorbar
colorbar_x = xleft + wspace + width + 0.01 
colorbar_y = ylower + (height-colorbar_height)/2 

# figure 1 data
Z=np.loadtxt("phase_matrix.txt", usecols=[2]).T
Z=Z.reshape(20,20)

p0=np.zeros((nx,ny))

# figure 2 data
fname = "p_0.dat"
p = np.loadtxt(fname,usecols=[2],unpack=True)
# enhanced by 16
p=p*60
pmin1 = min(p)
pmax1 = max(p)
p = p.reshape((nx,ny))
p = np.transpose(p)


# figure 3 data
fname2 ="p_3000_focusing.dat"
p2 = np.loadtxt(fname2,usecols=[2],unpack=True)
pmin2 = np.min(p2)
pmax2 = np.max(p2)
p2 = p2.reshape((nx,ny))
p2 = np.transpose(p2)

# figure 4 data
fname3 ="p_3000_no_focusing.dat"
p3 = np.loadtxt(fname3,usecols=[2],unpack=True)
# enhanced by 16
p3=p3*5
pmin3 = np.min(p3)
pmax3 = np.max(p3)
p3 = p3.reshape((nx,ny))
p3 = np.transpose(p3)


#figure 1 plot


ax1=fig.add_axes([xleft+0.015, ylower+hspace+0.025, width-0.015, height-0.04])
plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')

im=plt.pcolor(Z, edgecolors='k', linewidths=1,cmap='hsv')
#xticks([0,1,2,3,4,5,6],['1','20','40','60','80','90','101'])
#yticks([0,1,2,3,4,5,6],['1','20','40','60','80','90','101'])

#cbar=plt.colorbar(im, ticks=np.linspace(0,2*np.pi,5))
#cbar.ax.set_yticklabels([r'0',r'$\pi/2$',\
                     #r'$\pi$',r'$3\pi/2$',\
                     #r'$2\pi$'],fontsize=13,fontproperties=tnr)
#cbar.update_ticks() 
cbar=plt.colorbar(im, shrink=0.9)
cbar.set_ticks(np.linspace(0,2*np.pi,5))
cbar.set_ticklabels([r'$0$',r'$\frac{\pi}{2}$',\
                     r'$\pi$',r'$\frac{3\pi}{2}$',\
                     r'$2\pi$'])                    
plt.clim(0,2*np.pi)

#plt.xlabel('x',labelpad=5,fontsize=20,fontproperties=tnri)
#plt.ylabel('y',labelpad=5,fontsize=20,fontproperties=tnri)

plt.xticks([0,5,10,15,20],['1','25','50','75','101'],fontsize=13,fontproperties=tnr)
plt.yticks([0,5,10,15,20],['1','25','50','75','101'],fontsize=13,fontproperties=tnr)
#plt.xticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
#plt.yticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
plt.text(0.15,1.07,'(a)', color='black',fontsize=16, transform=ax1.transAxes, \
        fontproperties=tnr, ha='center',va='center')
#add label to the colorbar
plt.text(1.15,1.05,'Phase', color='black',fontsize=15, transform=ax1.transAxes, \
        fontproperties=tnr, ha='center',va='center')




#figure 2
ax2=fig.add_axes([xleft+wspace, ylower+hspace, width, height])
plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')
im1=plt.imshow(p, origin='lower', extent=[1,nx,1,ny],cmap='PuRd')
#plt.xlabel('x',labelpad=5,fontsize=20,fontproperties=tnri)
#plt.ylabel('y',labelpad=5,fontsize=20,fontproperties=tnri)
plt.xlim(1,nx)
plt.ylim(1,ny)
plt.clim(pmin2,pmax2)
plt.xticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
plt.yticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
#plt.xticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
#plt.yticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
plt.text(0.15,0.85,'(b)', color='black',fontsize=16, transform=ax2.transAxes, \
        fontproperties=tnr, ha='center',va='center')
      
plt.text(0.7,0.7,r'$\times$ 60', color='black',fontsize=16, transform=ax2.transAxes, \
        fontproperties=tnr, ha='center',va='center')         

# colorbar
ax = fig.add_axes([colorbar_x, ylower, colorbar_width, height+hspace-0.06])
plt.tick_params(direction='out', bottom='off', top='on', left='off', right='off')
ticks_min = pmin2 #np.around(pmin2,decimals=5)
ticks_max = pmax2 #np.around(pmax2,decimals=5)
ticks_range = np.linspace(ticks_min,ticks_max,5)
print ticks_range
plt.clim(ticks_min,ticks_max)
cbar = plt.colorbar( cax=ax, ticks=ticks_range)

cl=cbar.ax.set_yticklabels([r'0',r'8.2$\times$$10^{-3}$',\
                     r'1.6$\times$$10^{-2}$',r'2.5$\times$$10^{-2}$',\
                     r'3.3$\times$$10^{-2}$'],fontsize=13,fontproperties=tnr)
# add a label for the whole colorbar
plt.text(1.9,1.07,'            Probability', color='black',fontsize=16, transform=ax.transAxes, \
        fontproperties=tnr, ha='center',va='center')





#figure 3
ax3=fig.add_axes([xleft, ylower, width, height])
plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')
im2=plt.imshow(p2, origin='lower', extent=[1,nx,1,ny],cmap='PuRd')
plt.clim(pmin2,pmax2)
#plt.xlabel('x',labelpad=5,fontsize=20,fontproperties=tnri)
#plt.ylabel('y',labelpad=5,fontsize=20,fontproperties=tnri)
plt.xlabel(r'Site index ($x$-axis)',labelpad=7,fontsize=16,fontproperties=tnr)
plt.ylabel(r'Site index ($y$-axis)',labelpad=2,fontsize=16,fontproperties=tnr)
plt.xlim(1,nx)
plt.ylim(1,ny)
#plt.xticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
#plt.yticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
plt.xticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
plt.yticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
plt.text(0.15,0.85,'(c)', color='black',fontsize=16, transform=ax3.transAxes, \
        fontproperties=tnr, ha='center',va='center')



#figure 4
ax4=fig.add_axes([xleft+wspace, ylower, width, height])
plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')
im3=plt.imshow(p3, origin='lower', extent=[1,nx,1,ny],cmap='PuRd')
plt.clim(pmin2,pmax2)
#plt.xlabel('x',labelpad=7,fontsize=20,fontproperties=tnri)
#plt.ylabel('y',labelpad=2,fontsize=20,fontproperties=tnri)
plt.xlim(1,nx)
plt.ylim(1,ny)
plt.xticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
plt.yticks([1,25,50,75,101],['','','','',''],fontsize=14,fontproperties=tnr)
#plt.xticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
#plt.yticks([1,25,50,75,101],fontsize=14,fontproperties=tnr)
plt.text(0.15,0.85,'(d)', color='black',fontsize=16, transform=ax4.transAxes, \
        fontproperties=tnr, ha='center',va='center')
plt.text(0.7,0.7,r'$\times$ 5', color='black',fontsize=16, transform=ax4.transAxes, \
        fontproperties=tnr, ha='center',va='center') 






plt.savefig('t-matrix-focusing.png', bbox_inches='tight', dpi=300)
plt.savefig('t-matrix-focusing.pdf', bbox_inches='tight', dpi=300)
plt.show()
