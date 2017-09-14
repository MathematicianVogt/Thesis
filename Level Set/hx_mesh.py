import matplotlib.pyplot as plt
import math
from compare_exact_numerical import compare_exact_numerical_2D
from time_mesh import time_mesh
import numpy as np
from pde import *


class hx:
		# a<=x<=b
	# c<=y<=d
	#dx - spatial step in x
	#dy - spatial step in y
	#nx - number steps in x direction
	#nx - number steps in y direction
	def __init__(self,a,b,c,d,nx,ny,Tmax,nt,BCs,IC,phi,epsilon,mu):
		#stencil update i+.5,j+.5
		#time n+.5

		self.x = np.linspace(a,b,nx)
		self.y = np.linspace(c,d,ny)
		dyt= (self.y[1]-self.y[0])/2.0
		for i in range(0,len(self.y)):
			self.y[i]+=dyt
		self.grid_x = np.linspace(a,b,nx)
		self.grid_y=np.linspace(c,d,ny)
		self.x_list=self.x
		self.y_list=self.y
		self.xsize = len(self.x)
		self.ysize=len(self.y)
		self.xx,self.yy = np.meshgrid(self.x, self.y, indexing = 'ij')
		self.mesh = (self.xx,self.yy)
		self.time_mesh = time_mesh(Tmax,nt*2)
		self.dt = 2.0*(self.time_mesh.time_step())
		self.hx_sol=[]
		self.IC=IC
		self.add_ic()
		self.phi = phi_hx(phi,self)
		self.interface_grid = self.phi.set_irregular_regular_points()
		self.epsilon=epsilon
		self.mu=mu
		self.BC=BCs

	def get_interface(self):
		return (self.xx,self.yy,self.interface_grid)
	def plot_interface(self):
		plt.pcolor(self.xx, self.yy, self.interface_grid,cmap='cool')
		plt.show()

	
	def add_ic(self):
		x1=self.xsize
		x2=self.ysize
		IC_cond = np.zeros((x1,x2))
		for i in range(0,len(self.x)):
			for j in range(0,len(self.y)):	
				IC_cond[i,j] = self.IC(self.x_list[i],self.y_list[j])
		self.hx_sol.append(IC_cond)

	def enforce_boundary_conditons(self, t):
		#bc dictionary
		bc=self.BC
		top = bc["top"]
		left = bc["left"]
		right=bc["right"]
		bottom=bc["bottom"]

		new_sol_boundary_conditions_enforced=np.zeros((self.xsize,self.ysize))

		for i in range(0,self.xsize):
			for j in range(0,self.ysize):

				#left_bc
				if(i==0 and j>=0):
					new_sol_boundary_conditions_enforced[i,j] = left(self.y_list[j],t)
					
				#bottom BC
				if(j==0 and i>=0):
					new_sol_boundary_conditions_enforced[i,j] = bottom(self.x_list[i],t)

				#top BC
				if(j==self.ysize-1 and i>=0):
					new_sol_boundary_conditions_enforced[i,j] = top(self.x_list[i],t)
				#right bc
				if(i==self.xsize-1 and j>=0):
					new_sol_boundary_conditions_enforced[i,j] = right(self.y_list[j],t)

		return new_sol_boundary_conditions_enforced

	def x_list(self):
		return self.xx
	def y_list(self):
		return self.yy

	def one_d_x(self):
		return self.x
	def one_d_y(self):
		return self.y

	def return_point_on_mesh(self,i,j):
		
		if(i<0 or j <0):
			return (None,None)

		try:
			return (self.x[i],self.y[j])
		except:
			return (None,None)


	def size(self):
		return np.shape(self.mesh)
	#spatial steps
	def h(self):
		dx = (self.x[1]-self.x[0])
		dy=2.0*(self.y[1]-self.y[0])
		return (dx,dy)
	# def half_h(self):
	# 	dx = (self.x[1]-self.x[0])
	# 	dy=(self.y[1]-self.y[0])
	# 	return (dx,dy)
	def build_sol_regular(self,t,ez):
		(dx,dy) = self.h()
		mu=self.mu
		epsilon=self.epsilon
		dt =self.dt
		previous_hx = self.hx_sol[-1]
		hx = self.enforce_boundary_conditons(t)
		
		for i in range(1,len(self.x_list)-1):
			for j in range(0,len(self.y_list)-1):
				hx[i,j] = previous_hx[i,j] + (dt/(mu(self.x_list[i],self.y_list[j])*dy))*(ez[i,j+1] - ez[i,j])
		self.hx_sol.append(hx)

	def previous_sol(self):
		return self.hx_sol[-1]

class phi_hx:
	def __init__(self,phi,mesh):
		self.phi=phi
		self.mesh=mesh

	def set_irregular_regular_points(self):
		mesh=self.mesh
		a = mesh.size()
		print a
		xsize = a[1]
		ysize=a[2]
		mesh_numbering = np.zeros((xsize,ysize))

		a=mesh.size()
		xsize=a[1]
		ysize=a[2]
		for i in range(0,xsize):

			for j in range(0,ysize):
				mesh_numbering[i,j] = self.regular_or_irregular(i,j)
				#print (i,j)



		return mesh_numbering

	def evaluate(self,xi,yi):
		return self.phi(xi,yi)
	def evaluate_min(self,xi,yi):
		if(xi is None or yi is None ):
			return float("inf")
		else:
			return self.phi(xi,yi)
	def evaluate_max(self,xi,yi):
		if(xi is None  or yi is None):
			return -float("inf")
		else:
			return self.phi(xi,yi)

	def regular_or_irregular(self,i,j):
		phi_min = self.phi_min(i,j)
		phi_max = self.phi_max(i,j)


		# if j==0:
		# 	print (phi_min,phi_max,i,j)
		# 	print self.mesh.return_point_on_mesh(i,j)

		return int(phi_min*phi_max<= 0)

	def phi_min(self,i,j):
		

		p1=self.mesh.return_point_on_mesh(i,j+1)
		p2=self.mesh.return_point_on_mesh(i,j+2)
		p3=self.mesh.return_point_on_mesh(i,j)
		


		p1=self.evaluate_min(p1[0],p1[1])
		p2 = self.evaluate_min(p2[0],p2[1])
		p3=self.evaluate_min(p3[0],p3[1])

		return min(p1,p2,p3)


	def phi_max(self,i,j):
		
		p1=self.mesh.return_point_on_mesh(i,j+1)
		p2=self.mesh.return_point_on_mesh(i,j+2)
		p3=self.mesh.return_point_on_mesh(i,j)
		


		p1=self.evaluate_max(p1[0],p1[1])
		p2 = self.evaluate_max(p2[0],p2[1])
		p3=self.evaluate_max(p3[0],p3[1])

		return max(p1,p2,p3)
		

