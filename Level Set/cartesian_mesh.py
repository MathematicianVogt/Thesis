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
	
	def x(self):
		return self.xx
	def y(self):
		return self.yy

	def return_point_on_mesh(self,i,j):
		
		try:
			return (self.xx[i,j],self.yy[i,j])
		except:
			return (None,None)


	def size(self):
		return np.shape(self.mesh)
	def h(self):
		dx = 2.0*(self.x[1]-self.x[0])
		dy=2.0*(self.y[1]-self.y[0])
		return (dx,dy)
	def half_h(self):
		dx = (self.x[1]-self.x[0])
		dy=(self.y[1]-self.y[0])
		return (dx,dy)

x = cartesian_mesh(0.0,1.0,0.0,2.0,100,99)
print x.return_point_on_mesh(0,98)
print x.size()
