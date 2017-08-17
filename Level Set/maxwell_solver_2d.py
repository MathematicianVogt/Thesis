from cartesian_mesh import cartesian_mesh
from level_set_function import level_set_function
import matplotlib.pyplot as plt
import math
from compare_exact_numerical import compare_exact_numerical_2D
class maxwell_solver_2d_circle:
	def __init__(self,r,xo,yo,a,b,c,d,nx,ny):
		
		self.mesh  = cartesian_mesh(a,b,c,d,nx,ny)
		phi = lambda x,y:math.sqrt((x-xo)**2 +(y-yo)**2) - r
		self.phi = level_set_function(phi,self.mesh)
	def plot_mesh(self):
		mesh = self.mesh
		x = mesh.x()
		y=mesh.y()
		z = self.phi.set_irregular_regular_points(mesh)
		# print z
		# plt.pcolor(x, y, z,cmap='cool')
		# plt.show()
		a=compare_exact_numerical_2D(z,z,x,y)
		b=a.absolute_diff()
		a=compare_exact_numerical_2D(z,b,x,y)
		a.plot()


x = maxwell_solver_2d_circle(.5,0,0,-1,1,-1,1,100,100)
x.plot_mesh()