from ex_mesh import *
from ey_mesh import *
from ez_mesh import *
from hx_mesh import *
from hy_mesh import *
from hz_mesh import *


a=0
b=1
c=0
d=1
nx=500
ny=500
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


p1.plot_interface()
p2.plot_interface()
p3.plot_interface()
p4.plot_interface()
p5.plot_interface()
p6.plot_interface()

