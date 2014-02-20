import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')

ka, E = np.loadtxt("dispersion_90.dat").transpose()
plt.tick_params(direction='out', bottom='on', top='off', left='on', right='off')
plt.subplots_adjust(top=0.7,bottom=0.15)
plt.xlim(-np.pi, np.pi)
plt.xlabel(r'ka',fontsize=20,fontproperties=tnri)
#plt.ylabel(r'E(k) ', fontsize=20,fontproperties=tnri)
plt.text(-np.pi-0.9, -30, r'$E(k)$', fontsize=20,fontproperties=tnri, rotation='vertical')
plt.text(-np.pi-0.9, 37, '(arbitrary unit)', fontsize=20,fontproperties=tnr, rotation='vertical')

plt.xticks([-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi],\
      [r"", r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r''],\
      fontsize=20)
plt.text(-np.pi, -73, r'$-\pi$', fontsize=20, ha='center')
plt.text(np.pi, -73, r'$\pi$', fontsize=20, ha='center')
plt.ylim(-60, 70)
plt.yticks([-60, -40, -20, 0, 20, 40, 60], ['-60', '-40', '-20', '0', '20', '40', '60'], \
            fontsize=20, fontproperties=tnri)

plt.title('Dispersion curve of an exciton in 1D array', fontsize=20,fontproperties=tnr)
plt.plot(ka, E)
plt.savefig('dispersion_curve.pdf', bbox_inches='tight', dpi=300)
plt.savefig('dispersion_curve.png', bbox_inches='tight', dpi=300)
plt.show()
