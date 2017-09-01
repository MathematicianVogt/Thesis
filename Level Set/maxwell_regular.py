from cartesian_mesh import cartesian_mesh
from level_set_function import level_set_function
import matplotlib.pyplot as plt
import math
from compare_exact_numerical import compare_exact_numerical_2D
from time_mesh import time_mesh
import numpy as np
from pde import *


#assumption is that boundary conditions are direchlet 
class maxwell_reg:
	def __init__(self,a,b,c,d,nx,ny,Tmax,nt,BC_list,epsilon,mu,ic_list,source_tm,source_te):
		self.spatial_mesh = cartesian_mesh(a,b,c,d,nx,ny)
		self.time_mesh = time_mesh(Tmax,nt)
		self.time_mesh =self.time_mesh.mesh()
		self.dx,self.dy = self.spatial_mesh.h()
		self.dt = self.time_mesh[1]-self.time_mesh[0]
		self.BC=BC_list
		self.epsilon = epsilon
		self.mu=mu
		self.size = self.spatial_mesh.size()
		self.xsize=self.size[1]
		self.ysize=self.size[2]
		self.ic_list=ic_list
		self.x_list=self.spatial_mesh.one_d_x()
		self.y_list = self.spatial_mesh.one_d_y()
		self.source_tm = source_tm
		self.source_te= source_te
	#return the lambda function associated with the (i+1) equation i.e self.BC["top"][0] 
	#will return the top boundary condition for the first equation. 
	def get_boundary_condition(self,position,ith_equation):
		boundary_condition=self.BC[position]
		return boundary_condition[ith_equation]
	

	def initial_condition(self):
		ic_list=self.ic_list
		ic1=ic_list[0]
		ic2=ic_list[1]
		ic3=ic_list[2]

		initial_mesh_1 = np.zeros((self.xsize,self.ysize))
		initial_mesh_2 = np.zeros((self.xsize,self.ysize))
		initial_mesh_3 = np.zeros((self.xsize,self.ysize))

		x_list=self.spatial_mesh.one_d_x()
		y_list=self.spatial_mesh.one_d_y()
		for i in range(0,self.xsize):
			for j in range(0,self.ysize):
				initial_mesh_1[i,j] = ic1(x_list[i],y_list[j])
				initial_mesh_2[i,j] = ic2(x_list[i],y_list[j])
				initial_mesh_3[i,j] = ic3(x_list[i],y_list[j])

		return (initial_mesh_1,initial_mesh_2,initial_mesh_3)

	def yee_method_TM(self):
		hx=[]
		hy=[]
		ez=[]

		ic1,ic2,ic3 = self.initial_condition()

		(s1,s2,s3) = self.source_tm
		hx.append(ic1)
		hy.append(ic2)
		ez.append(ic3)

		dt=self.dt
		dx=self.dx
		dy=self.dy

		#half time step for first two, full time step for others 
		for t in range(1,len(self.time_mesh)):
			print "time:" + str(self.time_mesh[t])
			(sol1,sol2,sol3) = self.enforce_boundary_conditons(self.time_mesh[t])
			for i in range(0,self.xsize):
				for j in range(0,ysize):
					#first time step
					if t==1:
						try:
							sol1[i,j+1] = ic1[i,j+1] - (dt/(self.mu(self.x_list[i],self.y_list[j+1],self.time_mesh[t])*dy))*(ic3[i,j+2] - ic3[i,j]) + s1(self.x_list[i],self.y_list[j], self.time_mesh[t])
							
						except: 
							pass


						try:
							sol2[i+1,j] = ic2[i+1,j] + (dt/(self.mu(self.x_list[i+1],self.y_list[j],self.time_mesh[t])*dx))*(ic3[i+2,j] - ic3[i,j]) +  s2(self.x_list[i],self.y_list[j],self.time_mesh[t])




						except:
							pass
					
					else:
						hx_n=hx[-1]
						hy_n=hy[-1]
						ez_n=ez[-1]
						try:
							sol1[i,j+1] = hx_n[i,j+1] - (dt/(self.mu(self.x_list[i],self.y_list[j+1],self.time_mesh[t])*dy))*(ez_n[i,j+2] - ez_n[i,j]) +  s1(self.x_list[i],self.y_list[j],self.time_mesh[t])
						except:
							pass

						try:
							sol2[i+1,j] = hy_n[i+1,j] + (dt/(self.mu(self.x_list[i+1],self.y_list[j],self.time_mesh[t])*dx))*(ez_n[i+2,j] - ez_n[i,j]) + s2(self.x_list[i],self.y_list[j],self.time_mesh[t])
						except:
							pass

			for i in range(0,self.xsize):
				for j in range(0,ysize):
					
					if t==1:
						try:
							
							if(i-1<0 or j-1<0):
								raise Exception("Problem in X Y")
							else:
								sol3[i,j] = ic3[i,j] + (dt/(self.epsilon(self.x_list[i],self.y_list[j],self.time_mesh[t])*dx))*( ((ic2[i+1,j]  -ic2[i-1,j]))/dx   - ((ic1[i,j+1] -ic1[i,j-1])/dy)      ) + s3(self.x_list[i],self.y_list[j],self.time_mesh[t])
						

						except:
							pass
					else:
						hx_n=hx[-1]
						hy_n=hy[-1]
						ez_n=ez[-1]

						if(i-1<0 or j-1<0):
								raise Exception("Problem in X Y")
						else:

							sol3[i,j] = ez_n[i,j] + (dt/(self.epsilon(self.x_list[i],self.y_list[j],self.time_mesh[t])*dx))*( ((hy_n[i+1,j]  -hy_n[i-1,j]))/dx   - ((hx_n[i,j+1] -hx_n[i,j-1])/dy)      ) + s3(self.x_list[i],self.y_list[j],self.time_mesh[t])
						


			hx.append(sol1)
			hy.append(sol2)
			ez.append(sol3)
		return (hx,hy,ez)

	def enforce_boundary_conditons(self, t):
		bc=self.BC
		top = bc["top"]
		left = bc["left"]
		right=bc["right"]
		bottom=bc["bottom"]

		sol1=np.zeros((self.xsize,self.ysize))
		sol2=np.zeros((self.xsize,self.ysize))
		sol3=np.zeros((self.xsize,self.ysize))

		for i in range(0,self.xsize):
			for j in range(0,self.ysize):

				#left_bc
				if(i==0 and j>=0):
					sol1[i,j] = left[0](self.y_list[j],t)
					sol2[i,j] = left[1](self.y_list[j],t)
					sol3[i,j] = left[2](self.y_list[j],t)
				#bottom BC
				if(j==0 and i>=0):
					sol1[i,j] = bottom[0](self.x_list[i],t)
					sol2[i,j] = bottom[1](self.x_list[i],t)
					sol3[i,j] = bottom[2](self.x_list[i],t)
				#top BC
				if(j==self.ysize-1 and i>=0):
					sol1[i,j] = top[0](self.x_list[i],t)
					sol2[i,j] = top[1](self.x_list[i],t)
					sol3[i,j] = top[2](self.x_list[i],t)
				#right bc
				if(i==self.xsize-1 and j>=0):
					sol1[i,j] = right[0](self.y_list[j],t)
					sol2[i,j] = right[1](self.y_list[j],t)
					sol3[i,j] = right[2](self.y_list[j],t)

		return (sol1,sol2,sol3)


	def yee_method_TE(self):
		ex=[]
		ey=[]
		hz=[]

		ic1,ic2,ic3 = self.initial_condition()
		(s1,s2,s3) = self.source_tm
		ex.append(ic1)
		ey.append(ic2)
		hz.append(ic3)

		dt=self.dt
		dx=self.dx
		dy=self.dy

		#half time step for first two, full time step for others 
		for t in range(1,len(self.time_mesh)):
			print "time:" + str(self.time_mesh[t])
			(sol1,sol2,sol3) = self.enforce_boundary_conditons(self.time_mesh[t])
			for i in range(0,self.xsize):
				for j in range(0,self.ysize):
					#first time step
					if t==1:
						try:

							
							if (j-1<0):
								raise Exception("Out Of Grid y ")
							
							else:
								sol1[i+1,j] = ic1[i+1,j] + (dt/(self.epsilon(self.x_list[i+1],self.y_list[j],self.time_mesh[t])*dy))*(ic3[i+1,j+1] - ic3[i+1,j-1]) + s1(self.x_list[i],self.y_list[j])
							
						except:
							pass
							
						try:
							if(i-1<0):
								raise Exception(" Out OF grid X")
							else:
								sol2[i,j+1] = ic2[i,j+1] - (dt/(self.epsilon(self.x_list[i],self.y_list[j+1],self.time_mesh[t])*dx))*(ic3[i+1,j+1] - ic3[i-1,j+1]) + s2(self.x_list[i],self.y_list[j])
							




						except:
							pass
					
					else:
						ex_n=ex[-1]
						ey_n=ey[-1]
						hz_n=hz[-1]
						try:
							if (j-1<0):
								raise Exception("Out Of Grid y ")
								
							else:
								sol1[i+1,j] = ex_n[i+1,j] +(dt/(self.epsilon(self.x_list[i+1],self.y_list[j],self.time_mesh[t])*dy))*(hz_n[i+1,j+1] - hz_n[i+1,j-1]) +  s1(self.x_list[i],self.y_list[j])
						
						except:
							pass	
							
						try:
							if(i-1 < 0):
								raise Exception("Out Of Grid x ")
							else:
								sol2[i,j+1] = hz_n[i,j+1] - (dt/(self.epsilon(self.x_list[i],self.y_list[j+1],self.time_mesh[t])*dx))*(hz_n[i+1,j+1] - hz_n[i-1,j+1]) + s2(self.x_list[i],self.y_list[j])

						except:
							pass
						

			for i in range(0,self.xsize):
				for j in range(0,self.ysize):
					
					if t==1:
						try:
							


							sol3[i+1,j+1] = ic3[i+1,j+1] + (dt/(self.mu(self.x_list[i+1],self.y_list[j+1],self.time_mesh[t])*dx))*( ((ic1[i+1,j+2]  -ic1[i+1,j]))/dy   - ((ic2[i+2,j+1] -ic2[i,j+1])/dx)      ) +  s3(self.x_list[i],self.y_list[j])
						except:
							pass
					else:
						ex_n=ex[-1]
						ey_n=ey[-1]
						hz_n=hz[-1]

						try:
							sol3[i+1,j+1] = hz_n[i+1,j+1] + (dt/(self.mu(self.x_list[i+1],self.y_list[j+1],self.time_mesh[t])*dx))*( ((ex_n[i+1,j+1]  -ex_n[i+1,j]))/dy  - ((ey_n[i+2,j+1] -ey_n[i,j+1])/dy)      ) + s3(self.x_list[i],self.y_list[j])
						except:
							pass


			ex.append(sol1)
			ey.append(sol2)
			hz.append(sol3)


		return (ex,ey,hz)
	def source_tm(self,xi,yj):

		return self.source_tm(xi,yj)
	def source_te(self,xi,yj):
		return self.source_te(xi,yj)


	def plotsol_te(self,sol):

		sol1=sol[0]
		sol2=sol[1]
		sol3=sol[2]
		makeMovieMax(self.x_list,self.y_list,sol1,self.time_mesh,"Ex","x","y")
		makeMovieMax(self.x_list,self.y_list,sol2,self.time_mesh,"Ey","x","y")
		makeMovieMax(self.x_list,self.y_list,sol3,self.time_mesh,"Ez","x","y")









a = lambda x,t: 0
b = lambda x,y: 0
bc={}
bc_list = [a,a,a]
bc["top"]=bc_list
bc["left"]=bc_list
bc["right"]=bc_list
bc["bottom"]=bc_list
epsilon = lambda x,y,t : 1
mu = lambda x,y,t : 1
ic = [b,b,b ]
source_te=[lambda x,y,t:0,lambda x,y,t:0,lambda x,y,t:0]
source_tm =[lambda x,y,t:0,lambda x,y,t:0,lambda x,y,t:0]

x=maxwell_reg(0.0,1.0,0.0,1.0,100,100,1.0,10,bc,epsilon,mu,ic,source_tm,source_te)
sol =x.yee_method_TE()


#x.plotsol_te(sol)