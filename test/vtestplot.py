from numpy import *
from matplotlib.pyplot import *
import asciidata
from cometsuite import _AU

figure(1)
beta, vx, vy, vz = [array(x) for x in asciidata.open('vtest-1.dat')]
vej = sqrt(vx**2 + vy**2 + vz**2)
loglog(beta, vej, '.', hold=0, label='data')
b = 10**arange(-5.0, -0.9, 0.1)
plot(b, sqrt(b/4.08642), label='ideal')
setp(gca(), xlabel='Beta', ylabel='v (km/s)')
legend(loc='upper left')
show()

figure(2)
beta, vx, vy, vz, rx, ry, rz = [array(x) for x in asciidata.open('vtest-2.dat')]
vej = sqrt(vx**2 + vy**2 + vz**2)
rh = sqrt(rx**2 + ry**2 + rz**2)
v = array(zip(vx, vy, vz))
r = array(zip(rx, ry, rz))
costh = (v * r).sum(1) / vej / -rh
loglog(beta, vej, '.', hold=0, label='data')
loglog(beta, vej / costh, '.', label='corrected')
b = 10**arange(-5.0, -0.9, 0.1)
plot(b, sqrt(b/4.08642), label='v(th_sun=0)')
setp(gca(), xlabel='Beta', ylabel='v (km/s)')
legend(loc='upper left')
show()

figure(3)
beta, vx, vy, vz, rx, ry, rz = [array(x) for x in asciidata.open('vtest-3.dat')]
vej = sqrt(vx**2 + vy**2 + vz**2)
rh = sqrt(rx**2 + ry**2 + rz**2)
v = array(zip(vx, vy, vz))
r = array(zip(rx, ry, rz))
costh = (v * r).sum(1) / vej / -rh
loglog(beta, vej, '.', hold=0, label='data')
loglog(beta, vej / costh * sqrt(rh / _AU), '.', label='corrected')
b = 10**arange(-5.0, -0.9, 0.1)
plot(b, 0.5 * sqrt(b), label='v(th_sun=0, rh=1)')
setp(gca(), xlabel='Beta', ylabel='v (km/s)')
legend(loc='upper left')
show()


