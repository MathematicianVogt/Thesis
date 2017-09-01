import numpy as np
class time_mesh:
	def __init__(self,T,nt):
		self.t_mesh = np.linspace(0,T,nt)
		self.dt = self.t_mesh[1]-self.t_mesh[0]

	def time_step(self):
		return self.dt
	def half_time_step(self):
		return self.dt/2.0

	def get_time_location(self,i):
		return self.t_mesh[i]
	def size(self):
		return length(self.t_mesh)
	def mesh(self):
		return self.t_mesh
		

