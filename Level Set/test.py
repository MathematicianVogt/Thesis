from ex_mesh import *
from ey_mesh import *
from ez_mesh import *
from hx_mesh import *
from hy_mesh import *
from hz_mesh import *
import pylab as plt

a=0
b=1
c=0
d=1
nx=100
ny=100
Tmax=10
nt=100
BC = {}
basic_bc = lambda x,y:0
BC["top"] = basic_bc
BC["bottom"] = basic_bc 
BC["left"] = basic_bc
BC["right"] = basic_bc
IC = lambda x,y: 0 
xo = .5
yo=.5
r=.25
phi = lambda x,y : (x-xo)**2 + (y-yo)**2. - r**2


p1 = ex(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi)
p2 = ey(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi)
p3 = ez(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi)
p4 = hx(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi)
p5 = hy(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi)
p6 = hz(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi)


# p1.plot_interface()
# p2.plot_interface()
# p3.plot_interface()
# p4.plot_interface()
# p5.plot_interface()
# p6.plot_interface()

c1=p1.get_interface();
c2=p2.get_interface();
c3=p3.get_interface();
c4=p4.get_interface();
c5=p5.get_interface();
c6=p6.get_interface();


plt.figure(1)
plt.subplot(231)
plt.pcolor(c1[0], c1[1], c1[2],cmap='cool')
plt.title("E1")
plt.subplot(232)
plt.pcolor(c2[0], c2[1], c2[2],cmap='cool')
plt.title("E2")
plt.subplot(233)
plt.pcolor(c3[0], c3[1], c3[2],cmap='cool')
plt.title("E3")
plt.subplot(234)
plt.pcolor(c4[0], c4[1], c4[2],cmap='cool')
plt.title("B1")
plt.subplot(235)
plt.pcolor(c5[0], c5[1], c5[2],cmap='cool')
plt.title("B2")
plt.subplot(236)
plt.pcolor(c6[0], c6[1], c6[2],cmap='cool')
plt.title("B3")
plt.show()
