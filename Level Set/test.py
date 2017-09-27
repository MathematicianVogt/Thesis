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
<<<<<<< Updated upstream
Tmax=1
nt=2
=======
Tmax=10
nt=100000
>>>>>>> Stashed changes
BC = {}
basic_bc = lambda x,y:math.sin(math.pi*y)
BC["top"] = basic_bc
BC["bottom"] = basic_bc 
BC["left"] = basic_bc
BC["right"] = basic_bc
IC = lambda x,y: x**2 + y**2
xo = .5
yo=.5
r=.25
epsilon = lambda x,y : 1.0
mu = lambda x,y : 1.0
#phi = lambda x,y : (x-xo)**2 + (y-yo)**2. - r**2


phi = lambda x,y : x-.5
#phi = lambda x,y : y-.5


p1 = ex(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)
p2 = ey(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)
p3 = ez(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)
p4 = hx(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)
p5 = hy(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)
p6 = hz(a,b,c,d,nx,ny,Tmax,nt,BC,IC,phi,epsilon,mu)

t1 = time_mesh(Tmax,nt)
t2 =time_mesh(Tmax,2*nt)

a =solution(p1,p2,p3,p4,p5,p6,t1,t2)
a.solve_TE()
a.solve_TM()
maxwellplot=maxwell_plot(p1,p2,p3,p4,p5,p6,Tmax,nt)
maxwellplot.plot_contour("movielol")

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
