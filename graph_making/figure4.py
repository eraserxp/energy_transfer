from pylab import *
import matplotlib.font_manager as fm
from FortranFile import *


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')
figure()

subplots_adjust(hspace=0.47)

ax1=subplot(311,adjustable='box', aspect=0.00833)
tick_params(direction='out',bottom='on', top='off', left='on', right='off')
ka, energy0 = loadtxt('dispersion_0.dat').T
ka, energy90 = loadtxt('dispersion_90.dat').T
dim = ka.size
zeroenergy = zeros(dim)
xlim(-pi,pi)
ylim(-150,150)
ytick_array = np.linspace(-150,150,6)
tick_params(bottom='on', top='off', left='on', right='off')
xticks([-pi, -pi/2, 0, pi/2, pi], [r'-$\pi$', r'- $\frac{\pi}{2}$', r'0', r'$\frac{\pi}{2}$', r'$\pi$'],fontsize=13,fontproperties=tnr)
yticks(ytick_array,fontsize=13,fontproperties=tnr)


xlabel(r'k', labelpad=1, fontsize=16,fontproperties=tnri)
#ylabel(r'E(k) (kHz)',fontsize=17,fontproperties=tnr)
text(-4.30,-80,r'E(k)',fontsize=16,fontproperties=tnri,\
     ha='center',va='center',rotation='vertical')
text(-4.30,30,r'(kHz)',fontsize=16,fontproperties=tnr,\
     ha='center',va='center',rotation='vertical')     
#text()
text(-1.1, 46, r'$\theta=90^{\circ}$',color='red',horizontalalignment='center',va='center')
text(1.2, -78, r'$\theta=0^{\circ}$',color='blue',horizontalalignment='center',va='center')
text(0.0, -27, r'$\theta\approx 57.3^{\circ}$',color='green',horizontalalignment='center',va='center')
plot(ka, energy0, '--', ka, zeroenergy, ka, energy90, ':')
text(0.1, 0.87, '(a)', fontsize=16,fontproperties=tnr,transform=ax1.transAxes, horizontalalignment='center',verticalalignment='center')



ax2=subplot(312,adjustable='box', aspect=13.267)
tick_params(direction='out',bottom='on', top='off', left='on', right='off')
time, theta = loadtxt('flag6_theta_vs_time.txt').T
time=time*1.E6
tick_params(bottom='on', top='off', left='on', right='off')
xlim(0,3001)
ylim(0,90.2)
ylabel(r'$\theta$ (degree)',fontsize=16,fontproperties=tnr)

xticks([0, 600, 1200, 1800, 2400, 3000], fontsize=12,fontproperties=tnr)
yticks([ 0, 18, 36, 54, 72, 90 ],fontsize=12,fontproperties=tnr)
plot(time,theta, 'k', linewidth=2.0)
#text(1500, 75, '(b)', fontsize=17, horizontalalignment='center',verticalalignment='center')
text(0.1, 0.87, '(b)', fontsize=16,fontproperties=tnr,transform=ax2.transAxes, \
     horizontalalignment='center',verticalalignment='center')



ax3=subplot(313) #,adjustable='box', aspect=1.705)
tick_params(direction='out',bottom='on', top='off', left='on', right='off')

nrow=3001
ncolumn=1201
f = FortranFile('flag6.bin')
m1 = f.readRealArray(nrow, ncolumn, 'd').T

pmax=m1.max()
pmin=m1.min()

mpart=m1[100:700,:]

imshow(mpart,cmap='PuRd',origin='lower',aspect=1.9904)
clim(pmin,pmax*1.2)


xlabel(r'Time ($\mu$s)',fontsize=16,fontproperties=tnr)
ylabel('Site index',fontsize=16,fontproperties=tnr)
#text(1500, 130, '(c)', fontsize=17, horizontalalignment='center',verticalalignment='center',color='white')
xticks([0, 600, 1200, 1800, 2400, 3000],fontsize=12,fontproperties=tnr)
yticks([0, 120, 240, 360, 480, 600],['100','220','340','460','580','700'], fontsize=12, fontproperties=tnr)
text(0.1, 0.87, '(c)', fontsize=16,fontproperties=tnr,transform=ax3.transAxes, \
     horizontalalignment='center',verticalalignment='center')

savefig("figure4.png", bbox_inches='tight',dpi=300)
show()
