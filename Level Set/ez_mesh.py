import matplotlib.pyplot as plt
import math
from compare_exact_numerical import compare_exact_numerical_2D
from time_mesh import time_mesh
import numpy as np
from pde import *
import time


class ez:
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
		self.x_list=self.x
		self.y_list=self.y
		self.xsize = len(self.x)
		self.ysize=len(self.y)
		self.xx,self.yy = np.meshgrid(self.x, self.y, indexing = 'ij')
		self.mesh = (self.xx,self.yy)
		self.time_mesh = time_mesh(Tmax,nt)
		self.dt = self.time_mesh.time_step()
		self.ez_sol=[]
		self.IC=IC
		self.add_ic()
		self.phi = phi_ez(phi,self)
		self.interface_grid = self.phi.set_irregular_regular_points()
		self.epsilon=epsilon
		self.mu=mu
		self.BC=BCs
	def get_phi(self):
		return self.phi
	def get_sol(self):
		return (self.xx,self.yy,self.ez_sol)

	def get_interface(self):
		return (self.xx,self.yy,self.interface_grid)
	def plot_interface(self):
		plt.pcolor(self.xx, self.yy, self.interface_grid,cmap='cool')
		plt.show()

	
	def add_ic(self):
		x1=self.xsize
		x2=self.ysize
		IC_cond = np.zeros((x1,x2))
		for i in range(0,self.xsize):
			for j in range(0,self.ysize):	
				IC_cond[i,j] = self.IC(self.x_list[i],self.y_list[j])
		self.ez_sol.append(IC_cond)
		# print np.shape(IC_cond)
		# time.sleep(2)

	def enforce_boundary_conditons(self, t):
		#bc dictionary
		bc=self.BC
		top = bc["top"]
		left = bc["left"]
		right=bc["right"]
		bottom=bc["bottom"]
		# print (self.xsize,self.ysize)
		new_sol_boundary_conditions_enforced=np.zeros((self.xsize,self.ysize))


		for i in range(0,self.xsize):
			for j in range(0,self.ysize):
				# print self.xsize
				# print self.ysize
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
		#print np.shape(new_sol_boundary_conditions_enforced)
	
		return new_sol_boundary_conditions_enforced

	def x_list(self):
		return self.xx
	def y_list(self):
		return self.yy

	def one_d_x(self):
		return self.x
	def one_d_y(self):
		return self.y

	def return_point_on_mesh(self,i,j,half):
		
		if(i<0 or j <0):
			return (None,None)

		try:
			dx,dy = self.h()
			if half=="":
				return (self.x[i],self.y[j])
			elif half=="px":
				return (self.x[i] +dx/2.0  ,self.y[j])
			elif half=="py":
				return (self.x[i]  ,self.y[j] +dy/2.0)
			elif half=="nx":
				return (self.x[i] -dx/2.0  ,self.y[j])
			elif half=="ny":
				return (self.x[i]  ,self.y[j] -dy/2.0)
		except:
			return (None,None)


	def size(self):
		return np.shape(self.mesh)
	#spatial steps
	def h(self):
		dx = (self.x[1]-self.x[0])
		dy=(self.y[1]-self.y[0])
		return (dx,dy)
	# def half_h(self):
	# 	dx = (self.x[1]-self.x[0])
	# 	dy=(self.y[1]-self.y[0])
	# 	return (dx,dy)
	def previous_sol(self):
		return self.ez_sol[-1]

	def build_sol_regular(self,t,hx,hy):
		(dx,dy) = self.h()
		mu=self.mu
		epsilon=self.epsilon
		dt =self.dt
		previous_ez = self.ez_sol[-1]
		ez = self.enforce_boundary_conditons(t)
		for i in range(1,len(self.x_list)-1):
			for j in range(1,len(self.y_list)-1):
				ez[i,j] = previous_ez[i,j] + (dt/(epsilon(self.x_list[i],self.y_list[j])))*((hy[i,j] - hy[i-1,j])/dx  - (hx[i,j] - hx[i,j-1])/dy )
		self.ez_sol.append(ez)
class phi_ez:
	def __init__(self,phi,mesh):
		self.phi=phi
		self.mesh=mesh

	def set_irregular_regular_points(self):
		mesh=self.mesh
		a = mesh.size()
		#print a
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
		

		p1=self.mesh.return_point_on_mesh(i,j,"py")
		p2=self.mesh.return_point_on_mesh(i,j,"px")
		p3=self.mesh.return_point_on_mesh(i,j,"nx")
		p4=self.mesh.return_point_on_mesh(i,j,"ny")
		


		p1=self.evaluate_min(p1[0],p1[1])
		p2 = self.evaluate_min(p2[0],p2[1])
		p3=self.evaluate_min(p3[0],p3[1])
		p4=self.evaluate_min(p4[0],p4[1])

		return min(p1,p2,p3,p4)


	def phi_max(self,i,j):
		
		p1=self.mesh.return_point_on_mesh(i,j,"py")
		p2=self.mesh.return_point_on_mesh(i,j,"px")
		p3=self.mesh.return_point_on_mesh(i,j,"nx")
		p4=self.mesh.return_point_on_mesh(i,j,"ny")
		


		p1=self.evaluate_max(p1[0],p1[1])
		p2 = self.evaluate_max(p2[0],p2[1])
		p3=self.evaluate_max(p3[0],p3[1])
		p4=self.evaluate_max(p4[0],p4[1])

		return max(p1,p2,p3,p4)
	def gradient(self,i,j):
		x=self.x_list[i]
		y=self.y_list[j]
		(dx,dy) = self.mesh.h()
		phi_x = (self.phi(x+dx,y) + self.phi(x,y))/(2.0*dx)
		phi_y = (self.phi(x,y+dy) + self.phi(x,y))/(2.0*dy)

		return (phi_x,phi_y)
	def projection_and_get_angle(self,x1,y1,alpha):
		grad=self.gradient(x1,y1)
		x_star = (x1 + alpha*grad[0], y1 + alpha*grad[1])
		normal_at_interface = self.gradient(x_star[0],x_star[1])
		theta = math.acos((normal_at_interface[0])/math.sqrt(normal_at_interface[0]**2 + normal_at_interface[0]**2))
		return (theta,x_star)
	def theta_matrix(self,irregular_grid_point_matrix):

		theta_matrix = np.zeros(self.xsize,self.ysize)
		for i in range(0,self.xsize):
			for j in range(0,self.ysize):
				if (irregular_grid_point_matrix[i,j]==1):
					(theta,x_star) = self.projection_and_get_angle(self.x_list[i],self.y_list[j])
					theta_matrix[i,j] = theta
				else:
					pass

		return theta_matrix
	def curv(self,i,j):
		grad = self.gradient(i,j)
		(dx,dy) = self.mesh.h()
		phi_x=grad[0]
		phi_y=grad[1]
		phi_xx = (self.phi(self.x_list[i] - dx,self.y_list[j]) - 2.0*self.phi(self.x_list[i],self.y_list[j]) + self.phi(self.x_list[i] + dx,self.y_list[j]))/dx**2
		phi_yy = (self.phi(self.x_list[i],self.y_list[j]-dy) - 2.0*self.phi(self.x_list[i],self.y_list[j]) + self.phi(self.x_list[i] ,self.y_list[j] + dy))/dy**2
		phi_xy = (self.phi(self.x_list[i] + dx,self.y_list[j] +dy) -self.phi(self.x_list[i] + dx,self.y_list[j] -dy) -self.phi(self.x_list[i] - dx,self.y_list[j] +dy)  + self.phi(self.x_list[i] - dx,self.y_list[j] -dy)  )/(4.0*dx*dy)

		k = (phi_xx*phi_y**2 -2.0*phi_x*phi_y*phi_xy + phi_yy*phi_x**2)/(1e-15 +math.pow(phi_x**2 + phi_y**2),1.5)
		return k
	def normal(self,i,j):
		grad = self.gradient(i,j)
		phi_x = grad[0]
		phi_y=grad[1]
		norm = math.sqrt(1e-15 + phi_x**2 + phi_y**2)
		return (phi_x/norm,phi_y/norm)
		