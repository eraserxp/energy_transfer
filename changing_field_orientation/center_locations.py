from numpy import *
from pylab import *

x, y = loadtxt("center_locations.dat", usecols=[1,2], unpack=True)
#(a, b) = polyfit(x,y,1)
#angle = arctan(a)*180/pi
#angle = angle.round(decimals=1)
#yfit = polyval([a,b],x)
title("The motion of wavepacket in space")
plot(x,y, 'k.')
xlabel('x')
ylabel('y')
#plot(x,yfit,'r')
savefig("wavepacket_motion_kax_0.5pi_kay_0.5pi_theta_0_90_phi_0.pdf")
show()

