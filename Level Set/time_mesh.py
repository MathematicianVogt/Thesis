import numpy as np
class time_mesh:
	def __init__(self,T,nt):
		self.t_mesh = np.linespace(0,T,nt)
		self.dt = t_mesh[1]-t_mesh[0]

	def time_step(self):
		return 2.0*self.dt
	def half_time_step(self):
		return self.dt

	def get_time_location(self,i):
		return self.t_mesh[i]
	def size(self):
		return length(self.t_mesh)
		

