from cartesian_mesh import cartesian_mesh
from level_set_function import level_set_function
import matplotlib.pyplot as plt
import math
from compare_exact_numerical import compare_exact_numerical_2D
from time_mesh import time_mesh


#assumption is that boundary conditions are direchlet 
class maxwell_reg:
	def __init__(a,b,c,d,nx,ny,Tmax,nt,BC_list,epsilon,mu):
		self.spatial_mesh = cartesian_mesh(a,b,c,d,nx,ny)
		self.time_mesh = time_mesh(Tmax,nt)
		self.dx,self.dy = self.mesh.h()
		self.BC=BC_list
		self.epsilon = epsilon
		self.mu=mu
	#return the lambda function associated with the (i+1) equation i.e self.BC["top"][0] 
	#will return the top boundary condition for the first equation. 
	def get_boundary_condition(position,ith_equation):
		boundary_condition=self.BC[position]
		return boundary_condition[ith_equation]
	

	#Build A in AX=B
	def build_A(self,finite_dif_coeiff,node_locations):
		size=self.mesh.size()
		xsize=size[1]
		ysize=size[2]
		A=np.zeros((3*xsize,3*ysize))
	#Build B in AX=B
	def build_B(self,bounday_condtions):
		size=self.mesh.size()
		xsize=size[1]
		ysize=size[2]
		B=np.zeros((xsize,ysize))

