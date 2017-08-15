from pde import *
import math
import pylab as plt
import scipy.interpolate as inter
from transformcor import *

'''
Geneate Mesh and Set CFL condition for our problem.
'''
#set spatial domain (cartesian mesh)
a=0
b=1
c=0
d=1

#set time
t0=0.0
t1=10.0
nx=101
ny=101
n_time=100
spatial_step_x=(b-a)/float(n)
spatial_step_y=(d-c)/float(n)

#set material constants
e_plus = 1.0
e_minus=1.2
u_plus=1.0
u_minus=1.0

c_max = math.sqrt(1.0/min(e_plus,e_minus)*min(u_plus,u_minus))

e=1.0
u=1.0
dummy_dt = (t1-t0)/float(n_time)
#print dummy_dt

#our guess for time step to keep with CFL condition
#stab short for stability
stab = math.sqrt(spatial_step_x**2 + spatial_step_y**2)*(1.0/(dummy_dt*c_max))
#try to take biggest time step possible, while keeping CFL condition.
while(stab<1):
	n_time=n_time*2
	dummy_dt = (t1-t0)/float(n_time)
	stab = math.sqrt(spatial_step_x**2 + spatial_step_y**2)*(1.0/(dummy_dt*c_max))

spatial_step_t=dummy_dt

n_time = int((t1-t0)/spatial_step_x)+1

#construct spatial and time dimension.
x_list=onedinterval(a,b,nx)
y_list=onedinterval(c,d,ny)
t_list=onedinterval(t0,t1,n_time)

'''
Construct Our 2-D Interface
'''
#center of circle
xo=0
yo=0
#radius of circle
r=.5

#set the interval which will define our paramterterization
variable_paramter = onedinterval(0.0,1.0,1000)

x_func = lambda p : xo + r*math.cos(2.0*math.pi*p)
y_func = lambda p : yo + r*math.sin(2.0*math.pi*p)



def generate_paramterization(x_func,y_func,t_list):
	x_list=[]
	y_list=[]
	print 'info'
	print len(t_list)


	for i in t_list:
		x_list.append(x_func(i))
		y_list.append(y_func(i))


	print (len(x_list),len(y_list))
	return (x_list,y_list)

def generate_spline_2d(xval,yval):
	tck, u = inter.splprep([xval, yval], s=0)
	out = inter.splev(u, tck)
	return (out[0],out[1])


def generate_spline_2d_tangent(xval,yval):
	tck, u = inter.splprep([xval, yval], s=0)
	derv = inter.splev(u, tck,der=1)
	return (derv[0],derv[1])

def generate_spline_2d_normal(xval,yval):
	tck, u = inter.splprep([xval, yval], s=0)
	derv = inter.splev(u, tck,der=1)
	return (-derv[1],derv[0])

def generate_curveature(xval,yval):
	tck, u = inter.splprep([xval, yval], s=0)
	derv = inter.splev(u, tck,der=1)
	derv2 = inter.splev(u, tck,der=2)
	curv=[]
	for i in range(0,len(xval)):
		curv.append(math.fabs( (derv2[1][i]*derv[0][i] - derv2[0][i]*derv[1][i] ))/math.pow(derv[0][i]**2  +derv[1][i]**2,1.5) )
	return curv

#generate interface
(x_param,y_param) = generate_paramterization(x_func,y_func,variable_paramter)
#generate tangent vectors-indexed component wise by parameter variable
(t_x,t_y) = generate_spline_2d_tangent(x_param,y_param)
#generate normal vectors - indexed component wise by parameter variable
(n_x,n_y) = generate_spline_2d_normal(x_param,y_param)
#generate curvature at ever point  - indexed component wise by parameter variable
curvature = generate_curveature(x_param,y_param)
#plot interface
#plt.plot(x_param,y_param)
#plt.show()
theta = get_theta(t_x,t_y)





