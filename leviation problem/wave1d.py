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
#solving utt=uxx+s(u,x,t)
class wave1d:
	def __init__(self, a,b,Nx,Tmax,Nt,initial_postion,initial_velocity,left_bc,left_bc_type,right_bc,right_bc_type,source,initial_postion_func):
		self.a=a
		self.b=b
		self.initial_postion=initial_postion
		self.initial_velocity=initial_velocity
		self.left_bc=left_bc
		self.left_bc_type=left_bc_type
		self.right_bc=right_bc
		self.right_bc_type=right_bc_type
		(self.x,self.dx)=np.linspace(a,b,Nx,retstep=True)
		(self.t,self.dt)=np.linspace(0,Tmax,Nt,retstep=True)
		print self.dx
		print self.dt
		self.source=source
		self.u = None
		self.interp=None
	def d_source(self,u_in,x_i,t_n,dt):
		return dt*self.source(u_in,x_i,t_n)


		
	def buildSys(self,unplus1,un,unminus1,x,tn,first_step,velocity,c):
		sys=[]
		dx=self.dx
		dt=self.dt
		left_bc=self.left_bc
		right_bc=self.right_bc
		if(first_step):

			
			#left boundary condition
			if(self.left_bc_type=='D'):
				sys.append(unplus1[0] - self.left_bc(None,tn))
			elif(self.left_bc_type=='N'):
				sys.append((1.0/dt**2)*(unplus1[0] - 2*dt*velocity(self.x[0]) - 2*un[0] + unplus1[0]) -(c/dx**2)*(un[1]-2*dx*left_bc(un[0],tn) -2*unplus1[0]+ unplus1[1]) - source(un[0],self.x[0],tn) )

			#interior
			for i in range(1,len(self.x)-1):
				sys.append((1.0/dt**2)*(unplus1[i] - 2*dt*velocity(self.x[i]) - 2*un[i] + unplus1[i]) -(c/dx**2)*(unplus1[i-1] -2*unplus1[i]+ unplus1[i+1]) - source(un[i],self.x[i],tn) )
			#right boundary condition
			if(self.right_bc_type=='D'):
				sys.append(unplus1[len(self.x)-1] - self.right_bc(None,tn))
			elif(self.right_bc_type=='N'):
				N=len(self.x)-1
				sys.append((1.0/dt**2)*(unplus1[N] - 2*dt*velocity(self.x[N]) - 2*un[N] + unplus1[N]) -(c/dx**2)*(unplus1[N-1] -2*unplus1[N]+ 2*dx*right_bc(un[N],tn) + unplus1[N-1]) -source(un[N],self.x[N],tn))


		else:
			#left boundary condition
			if(self.left_bc_type=='D'):
				#print un[0]
				#print self.left_bc(None,tn)
				sys.append(unplus1[0] - self.left_bc(None,tn))
			elif(self.left_bc_type=='N'):
				sys.append((1.0/dt**2)*(unminus1[0]- 2*un[0] + unplus1[0]) -(c/dx**2)*(un[1]-2*dx*left_bc(un[0],tn) -2*unplus1[0]+ unplus1[1]) - source(un[0],self.x[0],tn) )

			#interior
			for i in range(1,len(self.x)-1):
				sys.append((1.0/dt**2)*(unminus1[i]  - 2*un[i] + unplus1[i]) -(c/dx**2)*(unplus1[i-1]-2*unplus1[i]+ unplus1[i+1]) - source(un[i],self.x[i],tn) )
			#right boundary condition
			if(self.right_bc_type=='D'):
				sys.append(unplus1[len(self.x)-1] - self.right_bc(None,tn))
			elif(self.right_bc_type=='N'):
				N=len(self.x)-1
				sys.append((1.0/dt**2)*(unminus1[N] - 2*un[N] + unplus1[N]) -(c/dx**2)*(unplus1[N-1] -2*unplus1[N]+ 2*dx*right_bc(un[N],tn) + unplus1[N-1]) -source(un[N],self.x[N],tn))

		return sys




	def build_func(self):
		f=lambda unplus1,un,unminus1,x,tn,first_step,velocity,c: self.buildSys(unplus1,un,unminus1,x,tn,first_step,velocity,c)
		return f

	def assembleSystemandSolve(self):
		
		if(initial_postion_func is True):
			self.initial_postion = np.vectorize(self.initial_postion)
			self.initial_velocity = np.vectorize(self.initial_velocity)
			int_sol=self.initial_postion(self.x)
			int_vel=self.initial_velocity(self.x)
		else:
			self.initial_velocity = np.vectorize(self.initial_velocity)
			int_sol=self.initial_postion
			int_vel=self.initial_velocity(self.x)

		list_ver= list(int_sol)
		list_vel=list(int_vel)
		u_sol=[list_ver]
		F=self.build_func()
		c=1.0
		#start building root finding system
		t=0
		for i in range(0,len(self.t)):
			
			sol=None
			if(i==0):
				uk=u_sol[-1]
				sol=optimize.root(F,uk,args=(uk,None,self.x,t,True,self.initial_velocity,c))
				sol=sol.x
			else:
				ukminus1=u_sol[-2]
				uk=u_sol[-1]
				sol=optimize.root(F,uk,args=(uk,ukminus1,self.x,t,False,self.initial_velocity,c))
				sol=sol.x




			t+=self.dt
			u_sol.append(list(sol))
		self.u=u_sol
		A=self.interpolat_func(u_sol)
		self.interp = interpolate.RectBivariateSpline(self.t,self.x,A)
		
		# print c
		# xx,tt= np.meshgrid(self.x,self.t)
		# fig = plt.figure()
		# ax = fig.gca(projection='3d')
		# surf = ax.plot_surface(xx, tt, A, cmap=cm.coolwarm,linewidth=0, antialiased=False)
		# plt.show()


		
		#self.makeMovie(self.x,u_sol,self.t,'Wave Equation Solution','x','t')
	def get_interp_timederv(t,x):
		return self.interp(t,x,1,0)
	def get_intrp_space_derv(t,x):
		return self.interp(t,x,0,1)

	def interpolat_func(self,usol):
		t=0
		x=self.a
		point_array=[]
		value_array=[]
		A=np.zeros((len(self.t),len(self.x)))
		for n in range(0,len(self.t)):
			x=self.a
			
			currentSol=usol[n]
			for i in range(0,len(self.x)):
				A[n,i]=currentSol[i]
		return A

	def get_solution(self):
		return self.u

	def makeMovie(self,xList, solList,tList,title,xlab,ylab):
		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=title, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=60, metadata=metadata)

		fig = plt.figure()
		l, = plt.plot([], [], 'k-o')

		

		with writer.saving(fig, title+ ".mp4", 100):
			for i in range(0,len(tList)):

				x0 = xList
				y0 =solList[i]

				#plt.xlim([-15,15])
				#plt.ylim([-2,2])
				plt.ylim(-2, 2)
				#plt.xlim(xList[0]-1,xList[len(xList)-1]+1)
				plt.xlim(0,1)
				plt.plot(x0,y0)
				plt.title(title + " Time = " + str(tList[i]))
				plt.xlabel(xlab)
				plt.ylabel(ylab)
				writer.grab_frame()
				plt.clf()	












# f=lambda x: math.sin(2*math.pi*x)
# g=lambda x: 0.0
# lbc=lambda u,t: 0.0
# rbc=lambda u,t: .5
# source = lambda u,x,t : 0.0
# wave=wave1d(0.0,1.0,100,1.0,1000,f,g,lbc,'D',rbc,'D',source)
# wave.assembleSystemandSolve()