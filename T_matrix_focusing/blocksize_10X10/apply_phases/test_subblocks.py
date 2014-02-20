import numpy as np
import matplotlib.pyplot as plt

TNX=101
TNY=101
bx = 20
by = 10

xticks=np.arange(1,TNX+bx, bx)
yticks=np.arange(1,TNY+by,by)
x, y = np.loadtxt('fort.46', usecols=(2,3)).T
plt.plot(x,y,'ro')
plt.xlim(1,TNX)
plt.ylim(TNY,1)
plt.xticks(xticks)
plt.yticks(yticks)
plt.grid(True)
plt.show()
