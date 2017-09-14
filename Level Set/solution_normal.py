class solution:
	def __init__(self,ex_mesh,ey_mesh,ez_mesh,hx_mesh,hy_mesh,hz_mesh,time_mesh,half_time_mesh):
		self.ex_mesh=ex_mesh
		self.ey_mesh=ey_mesh
		self.ez_mesh=ez_mesh
		self.hx_mesh=hx_mesh
		self.hy_mesh=hy_mesh
		self.hz_mesh=hz_mesh
		self.time_mesh=time_mesh
		self.half_time_mesh=half_time_mesh

	def solve_TE(self):
		ex = self.ex_mesh
		ey=self.ey_mesh
		hz = self.ez_mesh

		for n in range(1,self.time_mesh.size()):
			print "Time: " + str(self.time_mesh.get_time_location(n))
			pex = ex.previous_sol()
			pey= ey.previous_sol()
			phz = hz.previous_sol()
			t = self.time_mesh.get_time_location(n)
			ex.build_sol_regular(t,phz)
			ey.build_sol_regular(t,phz)
			hz.build_sol_regular(t,pex,pey)





	def solve_TM(self):
		hx=self.hx_mesh
		hy=self.hy_mesh
		ez=self.hz_mesh
		for n in range(1,self.time_mesh.size()):
			print "Time: " + str(self.time_mesh.get_time_location(n))
			phx = hx.previous_sol()
			phy= hy.previous_sol()
			pez = ez.previous_sol()
			t = self.time_mesh.get_time_location(n)
			hx.build_sol_regular(t,pez)
			hy.build_sol_regular(t,pez)
			ez.build_sol_regular(t,phx,phy)
