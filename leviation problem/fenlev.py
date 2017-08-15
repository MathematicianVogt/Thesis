from dolfin import *
import numpy as np
import math
class target(Expression):
	def eval(self, value, t):
	      value[0] = math.sin(2*math.pi*t[0])
	  def value_shape(self):
	      return (1,)
class delta(Expression):
	def eval(self, value, t):
	      value[0] = math.sin(2*math.pi*t[0])
	  def value_shape(self):
	      return (1,)
class lev:
	def __init__(r0,v0,Tmax,Nt,target_func,alpha,mass,charge):
		self.r0=r0
		self.v0=v0
		self.T=Tmax
		self.Nt=Nt
		self.target_func=target_func
		(self.t_list,self.dt)=np.linspace(0,Tmax,Nt,retstep=True)
		self.alpha=alpha
		self.mass=mass
		self.charge=charge

	def L(self,t,T):



	def target(self,t_list,f):
		target_arr=[]
		for i in range(0,len(t_list)):
			target_arr.append(f(target_arr[i]))
		return target_arr
	def solve_problem(self):
		timemesh=IntervalMesh(100, 0, self.T)
		spatialmesh = BoxMesh(-5.0, -5.0, -5.0, 5.0, 5.0,5.0, 50, 50, 50)
		V = VectorFunctionSpace(timemesh, "CG", 1,dim=3)
		W = VectorFunctionSpace(spatialmesh, "CG", 1, dim = 8)
		Z = MixedFunctionSpace([W,W,W,W,W,W,W,W,W,W,W])
		v=Function(Z)
		(phi,A1, A2,A3, p1,p2,p3,p4,m1,m2,m3) = split(v)
	

		
		bc_phi = DirichletBC(W.sub(0), 0.0, "on_boundary")
		bc_A1 = DirichletBC(W.sub(1), 0.0, "on_boundary")
		bc_A2 = DirichletBC(W.sub(2), 0.0, "on_boundary")
		bc_A3 = DirichletBC(W.sub(3), 0.0, "on_boundary")
		

		bcs = [bc_phi,bc_A1,bc_A2,bc_A3]


		t=0
		dt= self.T/1000.0

		while t<=self.T

			t+=dt











