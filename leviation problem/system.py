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
class system:
	def __init__(a,b,Tmax,nx,nt,c,inital_pos,initial_vel,x0,v0):
		self.a=a
		self.b=b
		self.Tmax=Tmax
		self.nx=nx
		self.nt=nt
		self.c=c
		self.x=np.linspace(self.a,self.b,self.nx)
		self.t=np.linspace(0.0,self.Tmax,nt)

	def build_func(self):
		f=lambda unplus1,un,unminus1,phiplus1,phin,phiminus1,vn,wn,x,tn,first_step,velocity,c,givenfunc: self.buildSys(unplus1,un,unminus1,phiplus1,phin,phiminus1,vn,wnx,tn,first_step,velocity,c,givenfuncic)
		return f

	def buildSys(self):







	def solveprob(self):
		phi=[]
		A=[]
		v=[]
		w=[]



		t=0
		for i in range(0,len(self.t)):
			
			sol=None
			if(i==0):
				phik=u_sol[-1]
				Ak=A[-1]
				vk=v[-1]
				wk=w[-1]
				sol=optimize.root(F,uk,args=(uk,None,self.x,t,True,self.initial_velocity,c))
				sol=sol.x
			else:
				phikminus1=phik[-2]
				phik=phik[-1]
				Akminus1=Ak[-2]
				Ak=Ak[-1]
				vk=v[-1]
				wk=w[-1]
				sol=optimize.root(F,uk,args=(uk,ukminus1,self.x,t,False,self.initial_velocity,c))
				sol=sol.x








			t+=self.dt