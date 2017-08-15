from dolfin import *
import numpy as np
import math
#solving laplace(u) = u_tt + f
#such that u=0 on boundary
#u(x,0) = initial condition
#u_t(x,0) = initial velocity 
# f force term
initial_velocity = Expression("sin(t)",t=0 ,degree =2)
f= Expression("0",t=0 ,degree =2)
initial_condition = Expression("sin(x[0]*x[1]*2*pi)",t=0 ,degree =2)
boundary_cond = Expression ("0", t=0,degree =2)
initial_condition.pi=math.pi


#1d mesh
#mesh=IntervalMesh(100, 0.0, 1.0)
#2d mesh
mesh = RectangleMesh(Point(0,0), Point(1, 1), 50, 50,)
#3d mesh
#mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0) ,20, 20, 20)

#define function space and zero boundary conditions
V = FunctionSpace(mesh, "CG", 1)
bc = DirichletBC(V , boundary_cond , " on_boundary ")

#define semi-discrete variables
un=project(initial_condition,V)
unm1=project(initial_condition,V)
unp1=TrialFunction( V )
temp=Function(V)
#define test function
v= TestFunction( V )


#final time and time step
T = 1.0
dt=.001
t = 0


#Set bilinear form for first time step
afirst = -(2.0/dt**2)*unp1*v*dx  
Lfirst = inner ( grad ( un ) , grad ( v ) )*dx -(2.0/dt)*initial_velocity*v*dx -(2.0/dt**2)*un*v*dx + f*v*dx


bc = DirichletBC(V , boundary_cond , " on_boundary ")
#set bilinear form for general problem(After first time step)
a=-(1.0/dt**2)*unp1*v*dx  
L=inner ( grad ( un ) , grad ( v ) )*dx +(1.0/dt**2)*unm1*v*dx -(2.0/dt**2)*un*v*dx + f*v*dx
Afirst = assemble( afirst )
A= assemble(a)

file = File("wavequation.pvd", "compressed")
#set up list of spatial solution at each time step, might be helpful in future code development not needed for now.
sol_list=[]
sol_list.append(initial_condition)
'''

ADD INTITIAL CONDITION TO PLOT.
AT THE MOMENT WE PLOT u(x,t1)....u(x,tn)
We would like u(x,t0)....u(x,tn) - have to cast the datastructure correctly


'''
#start marching in time
while t <= T+dt:
	print str(t)
	#we solve the problem in the first time step this way
	if t<=dt:
		b = assemble(Lfirst)
		boundary_cond.t = t
		f.t=t
		bc.apply(Afirst,b)
		solve(Afirst,temp.vector(),b)
		
		
		unm1.assign(un)
		un.assign(temp)
		sol_list.append(temp)
	#after time step we can just march in time.
	else:
		b = assemble(L)
		boundary_cond.t = t
		f.t=t
		bc.apply(A,b)
		solve(A,temp.vector(),b)
		unm1.assign(un)
		un.assign(temp)
		sol_list.append(temp)
	
	temp.rename('temp','temp')
	#save current time solution
	t += dt
	file << temp,t
	
 	



















