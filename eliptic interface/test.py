import pylab as plt
import numpy as np
from elliptic import *
import math


a=0
b=1
c=0
d=1
nx=100
ny=100


BC = {}
bc1 = lambda x: x**2 +1
bc2 = lambda x:x**2
bc3 = lambda x:x**2
bc4 = lambda x:x**2+1
BC["top"] = bc1
BC["bottom"] = bc2
BC["left"] = bc3
BC["right"] = bc4

xo = .5
yo=.5
r=.25
bplus=110.0
bminus=-51.0

phi = lambda x,y : math.sqrt((x-xo)**2 + (y-yo)**2) - r

#phi = lambda x,y : 1.0
def f(x,y):
	return 500*math.sin(x*y)

def v(x,y):
	return 1.0
def w(x,y):
	return 1.0

def sigma(x,y):
	return 0

def beta(x,y):
	bplus=10.0
	bminus=-5
	# if(x**2 + y**2 <= .25):
	# 	return bminus
	# else:
	# 	return bplus
	return 1.0




#phi = lambda x,y : x-.5
#phi = lambda x,y : y-.5


p1 = elliptic(a,b,c,d,nx,ny,BC,phi,beta,sigma,f,v,w)
#p1.build_sol_regular()
p1.plot_interface()

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
