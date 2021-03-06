import scipy
import numpy as np
#2d mesh
class cartesian_mesh:
	# a<=x<=b
	# c<=y<=d
	#dx - spatial step in x
	#dy - spatial step in y
	#nx - number steps in x direction
	#nx - number steps in y direction
	def __init__(self,a,b,c,d,nx,ny):
		self.x = np.linspace(a,b,nx)
		self.y = np.linspace(c,d,ny)
		self.xx,self.yy = np.meshgrid(self.x, self.y, indexing = 'ij')
		self.mesh = (self.xx,self.yy)
	
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
		dx = 2.0*(self.x[1]-self.x[0])
		dy=2.0*(self.y[1]-self.y[0])
		return (dx,dy)
	def half_h(self):
		dx = (self.x[1]-self.x[0])
		dy=(self.y[1]-self.y[0])
		return (dx,dy)


# x = cartesian_mesh(0.0,1.0,0.0,2.0,100,99)
# print x.return_point_on_mesh(0,98)
# print x.size()
