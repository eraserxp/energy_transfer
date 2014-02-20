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


#FIGURE A
ax1=plt.subplot(211)
plt.tick_params(direction='out', top='off', right='off',length=8,width=1.2,pad=8)
img2=plt.imread('k_umklapp_new.jpg')
ydim, xdim, colordim = np.shape(img2)
xtick_number = 6
xtick_position = np.linspace(0,xdim, xtick_number)
plt.imshow(img2,aspect='auto')
plt.xticks(xtick_position,['', '', '', '', '', ''],fontsize=16,fontproperties=tnr)
#plt.xticks([])
plt.xlim(0,xdim)

ytick_number=3
ytick_position = np.linspace(0,ydim, ytick_number)
plt.yticks(ytick_position, [r'$\pi$', '0', r'-$\pi$'],fontsize=28,fontproperties=tnr)
plt.ylim(0,ydim)
#plt.xlabel(r't ($ms$)',fontsize=17)
plt.ylabel('ka',fontsize=28,fontproperties=tnri)
ax1.text(0.9,0.8,'(a)',fontsize=24,fontproperties=tnr,color='black',transform=ax1.transAxes)



# FIGURE B
ax2=plt.subplot(212)
plt.tick_params(direction='out', top='off', right='off',length=8,width=1.2,pad=8)
img2=plt.imread('x_umklapp_new.jpg')
ydim, xdim, colordim = np.shape(img2)
xtick_number = 6
xtick_position = np.linspace(0,xdim, xtick_number)
plt.imshow(img2,aspect='auto')
plt.xticks(xtick_position,['0','400','800', '1200','1600','2000'],fontsize=23,fontproperties=tnr)
plt.xlim(0,xdim)
plt.xlabel(r'Time ($\mu$s)',fontsize=28,fontproperties=tnr)


plt.ylim(0,ydim)
ytick_number = 3
ytick_position = np.linspace(0,ydim, ytick_number)
plt.yticks(ytick_position,['600','500','400'],fontsize=23,fontproperties=tnr)
plt.ylabel('Site index',fontsize=28,fontproperties=tnr)
ax2.text(0.9,0.8,'(b)',fontsize=24,fontproperties=tnr, color='black',transform=ax2.transAxes)



plt.savefig('bloch_osci.pdf',bbox_inches='tight', dpi=300)
plt.show()
