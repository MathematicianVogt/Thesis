from ex_mesh import *
from ey_mesh import *
from ez_mesh import *
from hx_mesh import *
from hy_mesh import *
from hz_mesh import *
import pylab as plt
import numpy as np
from solution_normal import *
from time_mesh import *
from maxwell_plot import *

a=0
b=1
c=0
d=1
nx=100
ny=100
Tmax=1
nt=2
Tmax=10
nt=10
BC = {}
basic_bc = lambda x,y:1.0
BC["top"] = basic_bc
BC["bottom"] = basic_bc 
BC["left"] = basic_bc
BC["right"] = basic_bc
IC = lambda x,y: 5.0
xo = .5
yo=.5
r=.25
epsilon = lambda x,y : 1.0
mu = lambda x,y : 1.0
#phi = lambda x,y : (x-xo)**2 + (y-yo)**2. - r**2





BC_ex={}
BC_ey={}
BC_hz={}

IC_ex = lambda x,y: x+y
IC_ey = lambda x,y: x+y
IC_hz = lambda x,y: x+y


BC_ex["top"] = lambda x,t : x+1.0+t
BC_ex["bottom"] = lambda x,t : x+t
BC_ex["left"] = lambda y,t :y+t
BC_ex["right"] = lambda y,t :1.0+y+t

BC_ey["top"] = lambda x,t : x+1.0-t
BC_ey["bottom"] = lambda x,t : x-t
BC_ey["left"] = lambda y,t :y-t
BC_ey["right"] = lambda y,t :1.0+y-t


BC_hz["top"] = lambda x,t : x+1.0
BC_hz["bottom"] = lambda x,t : x
BC_hz["left"] = lambda y,t :y
BC_hz["right"] = lambda y,t :1.0+y








e=1.0
u=1.0
dummy_dt = (Tmax)/float(nt)
#print dummy_dt

spatial_step_x=(b-a)/float(nx)
spatial_step_y=(d-c)/float(ny)
spatial_step_t=None

#CHANGE LATER FOR CFL CONDITION TO BE HOLDING!!!!!!!!!!!!!!!!!!!!!!
c_max=1.0
#our guess for time step to keep with CFL condition
#stab short for stability
stab = math.sqrt(spatial_step_x**2 + spatial_step_y**2)*(1.0/(dummy_dt*c_max))
#try to take biggest time step possible, while keeping CFL condition.
while(stab<1):
	nt=nt*2
	dummy_dt = (Tmax)/float(nt)
	stab = math.sqrt(spatial_step_x**2 + spatial_step_y**2)*(1.0/(dummy_dt*c_max))

spatial_step_t=dummy_dt

nt = int((Tmax)/spatial_step_t)+100
print nt



phi = lambda x,y : x-.5
#phi = lambda x,y : y-.5


p1 = ex(a,b,c,d,nx,ny,Tmax,nt,BC_ex,IC_ex,phi,epsilon,mu)
p2 = ey(a,b,c,d,nx,ny,Tmax,nt,BC_ey,IC_ey,phi,epsilon,mu)
p3 = ez(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)
p4 = hx(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)
p5 = hy(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)
p6 = hz(a,b,c,d,nx,ny,Tmax,nt,BC_hz,IC_hz,phi,epsilon,mu)

t1 = time_mesh(Tmax,nt)
t2 =time_mesh(Tmax,2*nt)

a =solution(p1,p2,p3,p4,p5,p6,t1,t2)
a.solve_TE()
#a.solve_TM()
maxwellplot=maxwell_plot(p1,p2,p3,p4,p5,p6,Tmax,nt)
#maxwellplot.plot_contour("movielol")
maxwellplot.plot_contour_TE("Solution Movie")

# p1.plot_interface()
# p2.plot_interface()
# c2=p2.get_interface()
# print c2[2]
# p3.plot_interface()
# p4.plot_interface()
# p5.plot_interface()
# p6.plot_interface()

# c1=p1.get_interface();
# c2=p2.get_interface();
# c3=p3.get_interface();
# c4=p4.get_interface();
# c5=p5.get_interface();
# c6=p6.get_interface();


# # print np.count_nonzero(c1[2])
# # print np.count_nonzero(c2[2])
# # print np.count_nonzero(c3[2])
# # print np.count_nonzero(c4[2])
# # print np.count_nonzero(c5[2])
# # print np.count_nonzero(c6[2])



# # # plt.pcolor(c1[0],c1[1],c1[2],cmap='cool')
# # # plt.show()

#plot interfaces
# plt.figure(1)
# plt.subplot(231)
# plt.pcolor(c1[0], c1[1], c1[2],cmap='cool')
# plt.title("E1")
# plt.subplot(232)
# plt.pcolor(c2[0], c2[1], c2[2],cmap='cool')
# plt.title("E2")
# plt.subplot(233)
# plt.pcolor(c3[0], c3[1], c3[2],cmap='cool')
# plt.title("E3")
# plt.subplot(234)
# plt.pcolor(c4[0], c4[1], c4[2],cmap='cool')
# plt.title("H1")
# plt.subplot(235)
# plt.pcolor(c5[0], c5[1], c5[2],cmap='cool')
# plt.title("H2")
# plt.subplot(236)
# plt.pcolor(c6[0], c6[1], c6[2],cmap='cool')
# plt.title("H3")
# plt.show()
