import numpy as np
import math
import pylab as plt
from scipy import optimize
import time
import matplotlib.animation as manimation
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import interpolate
from scipy import integrate
from scipy import optimize

class lev:
	def __init__(x0,v0,a,b,Tmax,c,alpha,Nx,Nt):
		self.x0=x0,
		self.v0=v0
		self.a=a
		self.b=b
		self.Tmax=Tmax
		self.c=c
		self.alpha=alpha
		self.Nx=Nx
		self.Nt=Nt
		(self.x,self.dx)=np.linspace(a,b,Nx,retstep=True)
		self.t=np.linspace(0,Tmax,Nt)
		self.best_A=None
		self.best_Phi=None
		self.best_R=None
	def l2_norm_distance(position,target,time_interval):
		diff_array_squared=[]
		for i in range(0,len(position)):
			diff_array_squared.append((position[i]-target[i])**2)

		integral = integrate.simps(diff_array_squared,time_interval)
		return integral
	def l2_spatial(initial_A,x_interval):
		integral=integrate.simps(initial_A,x_interval)



	def target_list(f,time_list):
		f = np.vectorize(f)
		path=f(time_list)
		return path


	def ode_system(y,t,m,q,A,phi):
		
		dpdv=phi.get_intrp_space_derv(t,y)
		dadt=get_interp_timederv(t,y)
		m=1
		q=1
		dydx=[y[1],(1.0/m)*(-q*dpdv -q*dadt)]
		return dydx


	def cost(AIC,a,b,Tmax,nx,nt,x0,v0,a,b,Tmax,c,alpha):


		#set up A equation
		f=lambda x: math.sin(2*math.pi*x)
		p=lambda x: 10.0
		g=lambda x: 0.0
		lbc=lambda u,t: 0.0
		rbc=lambda u,t: 0.0
		source1 = lambda u,x,t : (math.pi/c)*
		A=wave1d(a,b,nx,1.0,nt,AIC,g,lbc,'D',rbc,'D',source,False)
		A.assembleSystemandSolve()

		phi=wave1d(a,b,nx,1.0,nt,p,g,lbc,'D',rbc,'D',source,True)
		phi.assembleSystemandSolve()




		sol = odeint(ode_system, [x0,v0], self.t, args=(A, phi))
		r=sol[:,0]

		self.best_A=A.get_solution()
		self.best_Phi=phi.get_solution()
		self.best_R=r

		return .5*l2_norm_distance(r,target_list(lambda t : t,self.t)) + (alpha/2.0)*l2_spatial(AIC,self.x)



	def solve_problem():
		a=self.a
		b=self.b
		c=self.c
		Tmax=self.Tmax
		nx=self.Nx
		nt=self.Nt
		x=self.x
		a0=lambda x: 1.0
		a0 = np.vectorize(a0)
		IC=a0(x)

		sol=optimize.minimize(cost,IC,(a,b,Tmax,nx,nt,x0,v0,a,b,Tmax,c,alpha))
		best_ic=sol.x
		np.savez('bigfile', best_A, best_Phi,best_R,best_ic)





problem=lev(0.0,0.0,0.0,1.0,1.0,1.0,.01,10,1000)

