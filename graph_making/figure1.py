import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')


# set up the fonts
#plt.rc('font', family='serif')
#plt.rc('font', serif='Times New Roman')


#compile and run the fortran code
#os.system('ifort -o plot -w -heap-arrays plots1.f90')
#os.system('/home/pxiang/run_jobs/wavepacket_1d/DC_plus_AC/make_plots/plot')




title = 'white background'

# make four images
k_figure="K_CPLong.ppm"
x_figure="X_CPLong.ppm"
tmp="tmp.ppm"

fig_A='k_momentum_kick_0_3_micron_second.jpg'
fig_B='k_momentum_kick_3_203_micron_second.jpg'
fig_C='x_momentum_kick_0_3_micron_second.jpg'
fig_D='x_momentum_kick_3_203_micron_second.jpg'

produce_A = " convert " + k_figure + " -crop 1001x201+0+0 " + tmp \
            + "; " + " convert " + tmp + " -quality 100 " + fig_A

produce_B = " convert " + k_figure + " -crop 1000x201+1001+0 " + tmp \
            + "; " + " convert " + tmp + " -quality 100 " + fig_B
            
produce_C = " convert " + x_figure + " -crop 1001x201+0+0 " + tmp \
            + "; " + " convert " + tmp + " -quality 100 " + fig_C

produce_D = " convert " + x_figure + " -crop 1000x201+1001+0 " + tmp \
            + "; " + " convert " + tmp + " -quality 100 " + fig_D
            
os.system(produce_A)
os.system(produce_B) 
os.system(produce_C)
os.system(produce_D) 

os.system(" rm " + tmp)          


   
plt.figure(figsize=(8.5,8.5))
gs=gridspec.GridSpec(2,2,width_ratios=[1,2])
gs.update(wspace=0.11,top=0.97,hspace=0.11,right=0.96)
# TITLE FOR THE WHOLE FIGURE
#plt.suptitle(title, fontsize=12)


#FIGURE A
ax1=plt.subplot(gs[0])

plt.tick_params(direction='out', top='off', right='off',length=8,width=1.2,pad=8)
img2=plt.imread(fig_A)
ydim, xdim, colordim = np.shape(img2)
xtick_number = 3
xtick_position = np.linspace(0,xdim, xtick_number)
plt.imshow(img2,aspect='auto')
#plt.xticks(xtick_position,['0', '1.5', '3'],fontsize=16)
#plt.xlim(0,xdim)
plt.xticks(xtick_position,[' ', ' ', ' '],fontsize=23,fontproperties=tnr)

ytick_number=5
ytick_position = np.linspace(0,ydim, ytick_number)
plt.yticks(ytick_position, [r'$\pi$', r'$\frac{\pi}{2}$', ' ', r'$\frac{\pi}{2}$', r'$\pi$'],fontsize=28,fontproperties=tnr)
ax1.text(-0.13,0.5,'0',fontsize=23,ha='center',va='center',transform=ax1.transAxes,fontproperties=tnr)
ax1.text(-0.23,0.755,'-',fontsize=23,ha='center',va='center',transform=ax1.transAxes,fontproperties=tnr)
ax1.text(-0.22,1.0,'-',fontsize=23,ha='center',va='center',transform=ax1.transAxes,fontproperties=tnr)
plt.ylim(0,ydim)
#plt.xlabel(r't ($ms$)',fontsize=17)
plt.text(-0.33,0.5, 'ka', ha='center', va='center', rotation='vertical', transform=ax1.transAxes, fontsize=28,fontproperties=tnri)
#plt.ylabel('k',position=(-0.17,0.5), ha='center', va='center', transform=ax1.transAxes, fontsize=32,fontproperties=tnri)
#ax1.text(0.7,0.85,'(a)',fontsize=18,color='black',transform=ax1.transAxes,fontproperties=tnr)

#get large tick marks
#for l in ax1.get_xticklines() + ax1.get_yticklines():
  #l.set_markersize(6)
  #l.set_markeredgewidth(1.2)


# FIGURE B
ax2=plt.subplot(gs[1])
plt.tick_params(direction='out', top='off', right='off',length=8,width=1.2,pad=8)
img2=plt.imread(fig_B)
ydim, xdim, colordim = np.shape(img2)
xtick_number = 3
xtick_position = np.linspace(0,xdim, xtick_number)
plt.imshow(img2,aspect='auto')
#plt.xticks(xtick_position,['3','103', '203'],fontsize=16)
#plt.xlim(0,xdim)
plt.xticks(xtick_position,[' ',' ', ' '],fontsize=23,fontproperties=tnr)
#get large tick marks
#for l in ax2.get_xticklines() + ax2.get_yticklines():
  #l.set_markersize(6)
  #l.set_markeredgewidth(1.2)

plt.ylim(0,ydim)
ytick_number=5
ytick_position = np.linspace(0,ydim, ytick_number)
plt.yticks(ytick_position, ['','','','',''],fontsize=23,fontproperties=tnr)
#plt.yticks([])
plt.ylabel('')
#ax2.text(0.9,0.85,'(b)',fontsize=17,color='black',transform=ax2.transAxes,fontproperties=tnr)

# insert phase map
ax_inside=plt.axes([0.747, 0.57, .2, .2])
image_inside=plt.imread('phasemap.png')
plt.imshow(image_inside,aspect='equal')
plt.setp(ax_inside, xticks=[], yticks=[])


#FIGURE C
ax3=plt.subplot(gs[2])
plt.tick_params(direction='out', top='off', right='off',length=8,width=1.2,pad=8)
img2=plt.imread(fig_C)
ydim, xdim, colordim = np.shape(img2)
xtick_number = 3
xtick_position = np.linspace(0,xdim, xtick_number)
plt.imshow(img2,aspect='auto')
plt.xticks(xtick_position,['0', '1.5', '3'],fontsize=23,fontproperties=tnr)
plt.xlim(0,xdim)
plt.text(0.8, -0.22, r'Time ($\mu$s)',fontsize=28,transform=ax3.transAxes,fontproperties=tnr)

ytick_number=5
ytick_position = np.linspace(0,ydim, ytick_number)
plt.yticks(ytick_position, [r'201', r'150', r'100', r'50', r'1'],fontsize=23,fontproperties=tnr)
plt.ylim(0,ydim)
#plt.xlabel(r't ($ms$)',fontsize=17)
plt.ylabel('Site index',fontsize=28,fontproperties=tnr)
#ax3.text(0.7,0.85,'(c)',fontsize=17,color='black',transform=ax3.transAxes,fontproperties=tnr)

#get large tick marks
#for l in ax3.get_xticklines() + ax3.get_yticklines():
  #l.set_markersize(6)
  #l.set_markeredgewidth(1.2)


# FIGURE D
ax4=plt.subplot(gs[3])
plt.tick_params(direction='out', top='off', right='off',length=8,width=1.2,pad=8)
img2=plt.imread(fig_D)
ydim, xdim, colordim = np.shape(img2)
xtick_number = 3
xtick_position = np.linspace(0,xdim, xtick_number)
plt.imshow(img2,aspect='auto')
plt.xticks(xtick_position,['3', '100', '200'],fontsize=23,fontproperties=tnr)
plt.xlim(0,xdim)


plt.ylim(0,ydim)
ytick_number=5
ytick_position = np.linspace(0,ydim, ytick_number)
plt.yticks(ytick_position, ['', '', '', '', ''],fontsize=23,fontproperties=tnr)
plt.ylabel('')

#get large tick marks
#for l in ax4.get_xticklines() + ax4.get_yticklines():
  #l.set_markersize(6)
  #l.set_markeredgewidth(1.2)
  
#ax4.text(0.9,0.85,'(d)',fontsize=17,color='black',transform=ax4.transAxes,fontproperties=tnr)

#plt.text(0.3,0.85,'Time', fontsize=17)

#savefile = "DC-" + DC_fluctuation + "%-" + "AC-" + laser_fluctuation + "%.pdf"


#savefile = 'black_background.pdf'
#savefile = 'figure1.pdf'
savefile = 'figure1.png'
plt.savefig(savefile,dpi=300)

plt.show()



#os.system("convert k_umklapp.ppm k_umklapp.jpg")
#os.system("convert x_umklapp.ppm x_umklapp.jpg")
#os.system( 'convert x_umklapp.ppm -crop 2001x201+0+400  tmp.ppm;' 
#            + 'convert  tmp.ppm -quality 100 x_umklapp.jpg')
#os.system('python2.7 figure4.py')
