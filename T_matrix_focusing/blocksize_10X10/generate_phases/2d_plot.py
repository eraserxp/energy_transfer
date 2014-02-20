from FortranFile import *
from pylab import *
while True:
  ntime=raw_input('Enter the time point:')
  if 0< int(ntime) <= 5000:
    filename='c_'+str(ntime)+'.bin'
    f=FortranFile(filename)
    nrow = 101
    ncolumn = 101
    x=f.readComplexArray(nrow, ncolumn)
    x_abs=abs(x)
    probability = x_abs**2
    #summation=sum(probability)
    imshow(probability,origin='lower',extent=[1,nrow,1,ncolumn])
    colorbar()
    show()
  else:
    print '************************************'
    print 'Outside the range of 1 to 5000'
    print 'Exit ... '
    break
