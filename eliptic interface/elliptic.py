import matplotlib.pyplot as plt
import math
import numpy as np
from pde import *
import time


class elliptic:
		# a<=x<=b
	# c<=y<=d
	#dx - spatial step in x
	#dy - spatial step in y
	#nx - number steps in x direction
	#nx - number steps in y direction
	def __init__(self,a,b,c,d,nx,ny,BCs,phi,beta,sigma,f,v,w):
		#stencil update i+.5,j+.5
		#time n+.5

		self.x = np.linspace(a,b,nx)
		self.y = np.linspace(c,d,ny)
		self.x_list=self.x
		self.y_list=self.y
		self.xx,self.yy = np.meshgrid(self.x, self.y, indexing = 'ij')
		self.xsize = len(self.x)
		self.ysize=len(self.y)
		self.mesh = (self.xx,self.yy)
		self.sol=[]
		self.phi = phi_e(phi,self)
		self.interface_grid = self.phi.set_irregular_regular_points()
		self.beta=beta
		self.BC=BCs
		self.sigma=sigma
		self.f=f
		self.v=v
		self.w=w

	def get_phi(self):
		return self.phi
	def get_sol(self):
		return (self.xx,self.yy,self.ex_sol)

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
		self.ex_sol.append(IC_cond)

	def enforce_boundary_conditons(self):
		#bc dictionary
		bc=self.BC
		top = bc["top"]
		left = bc["left"]
		right=bc["right"]
		bottom=bc["bottom"]

		new_sol_boundary_conditions_enforced=np.zeros((self.xsize,self.ysize))

		for i in range(0,self.xsize):
			for j in range(0,self.ysize):

				# #left_bc
				if(i==0 and j>=0):
					new_sol_boundary_conditions_enforced[i,j] = left(self.y_list[j])
					
				#bottom BC
				if(j==0 and i>=0):
					new_sol_boundary_conditions_enforced[i,j] = bottom(self.x_list[i])

				#top BC
				if(j==self.ysize-1 and i>=0):
					new_sol_boundary_conditions_enforced[i,j] = top(self.x_list[i])
				# #right bc
				if(i==self.xsize-1 and j>=0):
					new_sol_boundary_conditions_enforced[i,j] = right(self.y_list[j])

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
		dx = self.x[1]-self.x[0]
		dy=(self.y[1]-self.y[0])
		return (dx,dy)
	# def half_h(self):
	# 	dx = (self.x[1]-self.x[0])
	# 	dy=(self.y[1]-self.y[0])
	# 	return (dx,dy)
	def build_sol_regular(self):
		(dx,dy) = self.h()
		beta = self.beta
		interface_grid=self.interface_grid
		sol_with_boundary=self.enforce_boundary_conditons()
		A = np.zeros(((len(self.x_list)-1)**2,(len(self.y_list)-1)**2))
		b = np.zeros(((len(self.x_list)-1)**2,1))
		print A
		for i in range(1,len(self.x_list)-1):
			for j in range(1,len(self.y_list)-1):
				if (interface_grid[i,j]==1):
					self.irregular(A,b,i,j,dx,sol_with_boundary)
				else:
					self.regular(A,b,i,j,dx,sol_with_boundary)
		print A
		sol = np.linalg.solve(A, b)
		#print np.size(sol)


		major_sol = np.zeros((len(self.x_list)-2,len(self.y_list)-2))
		counter1=0
		counter2=0

		for i in range(0,(len(self.x_list)-1)**2):
			
			if counter2==len(self.x_list)-1:
				counter2=1
				counter1+=1	
			
			major_sol[:,:] = 2.0
			counter2+=1

		major_sol= np.flipud(major_sol)
		print major_sol
		print np.shape(major_sol)
		print np.shape(sol_with_boundary)
		sol_with_boundary[1:len(self.x_list)-1,1:len(self.y_list)-1] = major_sol
		print sol_with_boundary

	def regular(self,A,b,i,j,h,bound):
		k1 = self.k_transform(i,j,len(self.x_list))-1
		print k1
		k2 = self.k_transform(i-1,j,len(self.x_list))-1
		k3 = self.k_transform(i+1,j,len(self.x_list)) -1
		k4 = self.k_transform(i,j+1,len(self.x_list)) -1
		k5 = self.k_transform(i,j-1,len(self.x_list)) -1
		
		xi=self.x_list[i]
		yj=self.y_list[j]

		if(i==1 and j==1):

			A[k1,k1] = (1.0/h**2)*(-self.beta(xi + h/2.0,yj)-self.beta(xi,yj+h/2.0)) + self.sigma(xi,yj)
			A[k1,k2] =0.0
			A[k1,k3] =(1.0/h**2)*self.beta(xi + h/2.0,yj)
			A[k1,k4] = (1.0/h**2)*self.beta(xi ,yj+ h/2.0)
			A[k1,k5] =0.0

			b[k1,0] =self.f(xi,yj) - (1.0/h**2)*self.beta(xi-h/2.0 ,yj)*bound[i-1,j] - (1.0/h**2)*self.beta(xi ,yj- h/2.0)*bound[i,j-1]
			
		elif(i == len(self.x_list) and j==1):
			A[k1,k1] = (1.0/h**2)*(-self.beta(xi - h/2.0,yj)-self.beta(xi,yj+h/2.0)) + self.sigma(xi,yj)
			A[k1,k2] =(1.0/h**2)*self.beta(xi - h/2.0,yj)
			A[k1,k3] =0.0
			A[k1,k4] = (1.0/h**2)*self.beta(xi ,yj+ h/2.0)
			A[k1,k5] =0.0

			b[k1,0] =self.f(xi,yj) - (1.0/h**2)*self.beta(xi +h/2.0 ,yj)*bound[i+1,j] - (1.0/h**2)*self.beta(xi  ,yj - h/2.0)*bound[i,j-1]

		elif (i==1 and j==len(self.y_list)-2):
			A[k1,k1] = (1.0/h**2)*(-self.beta(xi + h/2.0,yj)-self.beta(xi,yj-h/2.0) ) + self.sigma(xi,yj)
			A[k1,k2] =0.0
			A[k1,k3] =(1.0/h**2)*self.beta(xi + h/2.0,yj)
			A[k1,k4] = 0.0
			A[k1,k5] =(1.0/h**2)*self.beta(xi ,yj- h/2.0)

			b[k1,0] =self.f(xi,yj) - (1.0/h**2)*self.beta(xi -h/2.0  ,yj)*bound[i-1,j] - (1.0/h**2)*self.beta(xi  ,yj +h/2.0)*bound[i,j+1]

		elif(i==len(self.x_list)-2, j==len(self.y_list-2)):
			A[k1,k1] = (1.0/h**2)*(self.beta(xi - h/2.0,yj)-self.beta(xi,yj-h/2.0) ) + self.sigma(xi,yj)
			A[k1,k2] =(1.0/h**2)*self.beta(xi - h/2.0,yj)
			A[k1,k3] =0.0
			A[k1,k4] = 0.0
			A[k1,k5] =(1.0/h**2)*self.beta(xi ,yj- h/2.0)
			b[k1,0] =self.f(xi,yj) - (1.0/h**2)*self.beta(xi +h/2.0  ,yj)*bound[i+1,j] - (1.0/h**2)*self.beta(xi  ,yj+h/2.0)*bound[i,j+1]
		elif (i==1):
			A[k1,k1] = (1.0/h**2)*(-self.beta(xi + h/2.0,yj)-self.beta(xi,yj+h/2.0)-self.beta(xi,yj-h/2.0) ) + self.sigma(xi,yj)
			A[k1,k2] =0.0
			A[k1,k3] =(1.0/h**2)*self.beta(xi + h/2.0,yj)
			A[k1,k4] = (1.0/h**2)*self.beta(xi ,yj+ h/2.0)
			A[k1,k5] =(1.0/h**2)*self.beta(xi ,yj- h/2.0)

			b[k1,0] =self.f(xi,yj)  - (1.0/h**2)*self.beta(xi -h/2.0  ,yj)*bound[i-1,j]
		elif (j==1):
			A[k1,k1] = (1.0/h**2)*(-self.beta(xi + h/2.0,yj)-self.beta(xi - h/2.0,yj)-self.beta(xi,yj+h/2.0) ) + self.sigma(xi,yj)
			A[k1,k2] =(1.0/h**2)*self.beta(xi - h/2.0,yj)
			A[k1,k3] =(1.0/h**2)*self.beta(xi + h/2.0,yj)
			A[k1,k4] = (1.0/h**2)*self.beta(xi ,yj+ h/2.0)
			A[k1,k5] =0.0
			b[k1,0] =self.f(xi,yj)  - (1.0/h**2)*self.beta(xi  ,yj-h/2.0)*bound[i,j-1]
		elif(i==len(self.x_list)-2):
			A[k1,k1] = (1.0/h**2)*(-self.beta(xi + h/2.0,yj)-self.beta(xi - h/2.0,yj)-self.beta(xi,yj-h/2.0) ) + self.sigma(xi,yj)
			A[k1,k2] =(1.0/h**2)*self.beta(xi - h/2.0,yj)
			A[k1,k3] =0.0
			A[k1,k4] = (1.0/h**2)*self.beta(xi ,yj+ h/2.0)
			A[k1,k5] =(1.0/h**2)*self.beta(xi ,yj- h/2.0)
			b[k1,0] =self.f(xi,yj)  - (1.0/h**2)*self.beta(xi +h/2.0  ,yj)*bound[i+1,j]
		elif(j==len(self.y_list)-2):
			A[k1,k1] = (1.0/h**2)*(-self.beta(xi + h/2.0,yj)-self.beta(xi - h/2.0,yj)-self.beta(xi,yj+h/2.0)-self.beta(xi,yj-h/2.0) ) + self.sigma(xi,yj)
			A[k1,k2] =(1.0/h**2)*self.beta(xi - h/2.0,yj)
			A[k1,k3] =(1.0/h**2)*self.beta(xi + h/2.0,yj)
			A[k1,k4] = (1.0/h**2)*self.beta(xi ,yj+ h/2.0)
			A[k1,k5] =(1.0/h**2)*self.beta(xi ,yj- h/2.0)
			b[k1,0] =self.f(xi,yj)  - (1.0/h**2)*self.beta(xi  ,yj+h/2.0)*bound[i,j+1]
		else:
			A[k1,k1] = (1.0/h**2)*(-self.beta(xi + h/2.0,yj)-self.beta(xi - h/2.0,yj)-self.beta(xi,yj+h/2.0)-self.beta(xi,yj-h/2.0) ) + self.sigma(xi,yj)
			A[k1,k2] =(1.0/h**2)*self.beta(xi - h/2.0,yj)
			A[k1,k3] =(1.0/h**2)*self.beta(xi + h/2.0,yj)
			A[k1,k4] = (1.0/h**2)*self.beta(xi ,yj+ h/2.0)
			A[k1,k5] =(1.0/h**2)*self.beta(xi ,yj- h/2.0)
			b[k1,0] =self.f(xi,yj)



	def irregular(self,A,b,i,j,h,bound):
		k1 = self.k_transform(i,j,len(self.x_list))
		k2 = self.k_transform(i-1,j,len(self.x_list)) 
		k3 = self.k_transform(i+1,j,len(self.x_list)) 
		k4 = self.k_transform(i,j+1,len(self.x_list)) 
		k5 = self.k_transform(i,j-1,len(self.x_list)) 

		# A[k1,k1] = 
		# A[k1,k2] =
		# A[k1,k3] =
		# A[k1,k4] = 
		# A[k1,k5] =
		# b[k1,0] = 

	def build_irregular_weights(self):
		pass


	def k_transform(self,j,i,real_length_x):
		return i+ (real_length_x-1)*(j-1)


	def previous_sol(self):
		return self.sol[-1]
class phi_e:
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
		
		(dx,dy)  = self.mesh.h()
		p1=self.mesh.return_point_on_mesh(i,j)
		p2=self.mesh.return_point_on_mesh(i,j+1)
		p3=self.mesh.return_point_on_mesh(i,j-1)
		p4=self.mesh.return_point_on_mesh(i+1,j)
		p5=self.mesh.return_point_on_mesh(i-1,j)
		#(theta,x_star ) = self.projection_and_get_angle(self.mesh.x_list[i],self.mesh.y_list[j],.001) 
		#p6=(x_star[0] + dx , x_star[1] + dy)
		p6=self.mesh.return_point_on_mesh(i+1,j+1)
		




		p1=self.evaluate_min(p1[0],p1[1])
		p2 = self.evaluate_min(p2[0],p2[1])
		p3=self.evaluate_min(p3[0],p3[1])
		p4=self.evaluate_min(p4[0],p4[1])
		p5=self.evaluate_min(p5[0],p5[1])
		p6=self.evaluate_min(p6[0],p6[1])
		return min(p1,p2,p3,p4,p5,p6)

	def phi_max(self,i,j):
		(dx,dy)  = self.mesh.h()
		p1=self.mesh.return_point_on_mesh(i,j)
		p2=self.mesh.return_point_on_mesh(i,j+1)
		p3=self.mesh.return_point_on_mesh(i,j-1)
		p4=self.mesh.return_point_on_mesh(i+1,j)
		p5=self.mesh.return_point_on_mesh(i-1,j)
		p6=self.mesh.return_point_on_mesh(i+1,j+1)
		#(theta,x_star ) = self.projection_and_get_angle(self.mesh.x_list[i],self.mesh.y_list[j],.001) 
		#p6=(x_star[0] + dx , x_star[1] + dy)



		p1=self.evaluate_max(p1[0],p1[1])
		p2 = self.evaluate_max(p2[0],p2[1])
		p3=self.evaluate_max(p3[0],p3[1])
		p4=self.evaluate_max(p4[0],p4[1])
		p5=self.evaluate_max(p5[0],p5[1])
		p6=self.evaluate_max(p6[0],p6[1])
		return max(p1,p2,p3,p4,p5,p6)

	def gradient(self,i,j):
		x=self.mesh.x_list[i]
		y=self.mesh.y_list[j]
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



