#from numpy import *
#from pylab import *
from pylab import *
import matplotlib.font_manager as fm


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')

# for the case where theta is changing from 0 to 90 degree and phi is fixed
def theta0(time):
  if 0<=time<=1.0:
    return 0.0
  elif (1.0<time<2.0):
    return 90.0*sin(pi*(time-1.0)/2)**2
  else:
    return 90.0
    
def phi0(time):
  return 0.0      

#for the case where phi is changing from 0 to 90 degree and theta is fixed
def phi1(time):
  if 0<=time<=1.0:
    return 0.0
  elif (1.0<time<2.0):
    return 90.0*sin(pi*(time-1.0)/2)**2
  else:
    return 90.0
    
def theta1(time):
  return 0.0  
  
timelist=linspace(0,3.0,1000)
theta0_list=[theta0(t) for t in timelist]
phi0_list=[phi0(t) for t in timelist]
theta1_list=[theta1(t) for t in timelist]
phi1_list=[phi1(t) for t in timelist]


figure()
subplots_adjust(hspace=0.35,wspace=0.45)
    
ax1=subplot(221)
tick_params(direction='out',bottom='on', top='off', left='on', right='off')
x, y = loadtxt("center_locations_theta_changing.txt", usecols=[1,2], unpack=True)
plot(x,y, 'k.')
xticks([-250,-125,0,125,250],fontsize=13,fontproperties=tnr)
yticks([-250,-125,0,125,250],fontsize=13,fontproperties=tnr)
xlabel(r'Site index ($x$-axis)',fontsize=16,fontproperties=tnr)
ylabel(r'Site index ($y$-axis)',fontsize=16,fontproperties=tnr)
xlim(-250,250)
ylim(-250,250)
text(0,30,r'$t=0$',fontsize=15,ha='center')
text(-7,-220,r'$t=3$ ms',fontsize=15,fontproperties=tnr)
text(0.1, 0.8, '(a)', fontsize=16,fontproperties=tnr,transform=ax1.transAxes) #, horizontalalignment='center',verticalalignment='center')


ax2=subplot(222)
tick_params(direction='out',bottom='on', top='off', left='on', right='off')
x, y = loadtxt("center_locations_phi_changing.txt", usecols=[1,2], unpack=True)
plot(x,y, 'k.')
xticks([-250,-125,0,125,250],fontsize=13,fontproperties=tnr)
yticks([-250,-125,0,125,250],fontsize=13,fontproperties=tnr)
xlabel(r'Site index ($x$-axis)',fontsize=16,fontproperties=tnr)
ylabel(r'Site index ($y$-axis)',fontsize=16,fontproperties=tnr)
xlim(-250,250)
ylim(-250,250)
text(0,30,r'$t=0$',fontsize=15,ha='center')
text(61,110,r'$t=3$ ms',fontsize=15,fontproperties=tnr)
text(0.1, 0.8, '(b)', fontsize=16,fontproperties=tnr,transform=ax2.transAxes)


ax3=subplot(223)
tick_params(direction='out',bottom='on', top='off', left='on', right='off')
plot(timelist, theta0_list, 'r', timelist, phi0_list, 'b')
ylim(-4,94)
xticks([0,0.5,1,1.5,2,2.5,3],['0','0.5','1.0','1.5','2.0','2.5','3.0'],fontsize=13,fontproperties=tnr)
yticks([0,15,30,45,60,75,90],fontsize=13,fontproperties=tnr)
xlabel(r'Time (ms)',fontsize=16,fontproperties=tnr)
ylabel(r'Angle (degree)',fontsize=16,fontproperties=tnr)
text(1.27,48,r'$\theta$',color='red',fontsize=16)
text(2,6.0,r'$\phi$',color='blue',fontsize=16)
text(0.1, 0.8, '(c)', fontsize=16,fontproperties=tnr,transform=ax3.transAxes)

ax4=subplot(224)
tick_params(direction='out',bottom='on', top='off', left='on', right='off')
plot(timelist, theta1_list, 'r', timelist, phi1_list, 'b')
ylim(-4,94)
xticks([0,0.5,1,1.5,2,2.5,3],['0','0.5','1.0','1.5','2.0','2.5','3.0'],fontsize=13,fontproperties=tnr)
yticks([0,15,30,45,60,75,90],fontsize=13,fontproperties=tnr)
xlabel(r'Time (ms)',fontsize=16,fontproperties=tnr)
ylabel(r'Angle (degree)',fontsize=16,fontproperties=tnr)
text(1.27,48,r'$\phi$',color='blue',fontsize=16)
text(2,6.0,r'$\theta$',color='red',fontsize=16)
text(0.1, 0.8, '(d)', fontsize=16,fontproperties=tnr,transform=ax4.transAxes)

savefig("control-2d-wavepacket.pdf", bbox_inches='tight',dpi=300)
show()

