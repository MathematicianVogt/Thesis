import numpy as np

class level_set_function:
	def __init__(self,phi_func,mesh):
		self.phi=phi_func
		self.mesh = mesh
	def set_irregular_regular_points(self,mesh):
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
		xmm = self.mesh.return_point_on_mesh(i-2,j)
		xm = self.mesh.return_point_on_mesh(i-1,j)
		xpp = self.mesh.return_point_on_mesh(i+2,j)
		xp = self.mesh.return_point_on_mesh(i+1,j)
		ymm=self.mesh.return_point_on_mesh(i,j-2)
		ym=self.mesh.return_point_on_mesh(i,j-1)
		ypp=self.mesh.return_point_on_mesh(i,j+2)
		yp=self.mesh.return_point_on_mesh(i,j+1)
		
		phixmm=self.evaluate_min(xmm[0],xmm[1])
		phixm = self.evaluate_min(xm[0],xm[1])
		phixpp=self.evaluate_min(xpp[0],xpp[1])
		phixp = self.evaluate_min(xp[0],xp[1])
		phiymm=self.evaluate_min(ymm[0],ymm[1])
		phiym = self.evaluate_min(ym[0],ym[1])
		phiypp=self.evaluate_min(ypp[0],ypp[1])
		phiyp = self.evaluate_min(yp[0],yp[1])
		return min(phixmm,phixm, phixpp,phixp,phiymm,phiym,phiypp,phiyp)

	def phi_max(self,i,j):
		xmm = self.mesh.return_point_on_mesh(i-2,j)
		xm = self.mesh.return_point_on_mesh(i-1,j)
		xpp = self.mesh.return_point_on_mesh(i+2,j)
		xp = self.mesh.return_point_on_mesh(i+1,j)
		ymm=self.mesh.return_point_on_mesh(i,j-2)
		ym=self.mesh.return_point_on_mesh(i,j-1)
		ypp=self.mesh.return_point_on_mesh(i,j+2)
		yp=self.mesh.return_point_on_mesh(i,j+1)
		
		phixmm=self.evaluate_max(xmm[0],xmm[1])
		phixm = self.evaluate_max(xm[0],xm[1])
		phixpp=self.evaluate_max(xpp[0],xpp[1])
		phixp = self.evaluate_max(xp[0],xp[1])
		phiymm=self.evaluate_max(ymm[0],ymm[1])
		phiym = self.evaluate_max(ym[0],ym[1])
		phiypp=self.evaluate_max(ypp[0],ypp[1])
		phiyp = self.evaluate_max(yp[0],yp[1])
		if j==0:
			print ymm
			print ym
			print [phixmm,phixm, phixpp,phixp,phiymm,phiym,phiypp,phiyp]
		return max(phixmm,phixm, phixpp,phixp,phiymm,phiym,phiypp,phiyp)


## taking perspective h'=2h
## then h'/2 = h. This way we can index by integers. 






