import pylab as plt
import numpy as np
class maxwell_plot:
	def __init__(self,ex,ey,ez,hx,hy,hz,Tmax,nt):
		self.ex=ex
		self.ey=ey
		self.ez=ez
		self.hx=hx
		self.hy=hy
		self.hz=hz
		t_list = np.linspace(0,Tmax,nt):
		t_list_shit = [x + (t_list[1]-t_list[0])/2.0 for x in t_list]
		

	def plot_contour(self):
		(x1,y1,ex_sol) = self.ex.get_sol()
		ex_interface_func = self.ex.get_phi()

		(x2,y2,ey_sol) = self.ey.get_sol()
		ey_interface_func = self.ey.get_phi()

		(x3,y3,ez_sol) = self.ez.get_sol()
		ez_interface_func = self.ez.get_phi()

		

		(x4,y4,hx_sol) = self.hx.get_sol()
		hx_interface_func = self.hx.get_phi()

		(x5,y5,hy_sol) = self.hy.get_sol()
		hy_interface_func = self.hy.get_phi()

		(x5,y5,hz_sol) = self.hz.get_sol()
		hz_interface_func = self.hz.get_phi()
		plt.figure()
		
		plt.subplot(231)
		plt.contour(x1,y1,ex_sol,cmap='cool')
		plt.contour(x1,y1,ex_interface_func(x1,y1))

		plt.subplot(232)
		plt.contour(x2,y2,ey_sol,cmap='cool')
		plt.contour(x2,y2,ey_interface_func(x2,y2))

		plt.subplot(233)
		plt.contour(x3,y3,ez_sol,cmap='cool')
		plt.contour(x3,y3,ez_interface_func(x3,y3))


		plt.subplot(234)
		plt.contour(x4,y4,hx_sol,cmap='cool')
		plt.contour(x4,y4,hx_interface_func(x4,y4))


		plt.subplot(235)
		plt.contour(x5,y5,hy_sol,cmap='cool')
		plt.contour(x5,y5,hy_interface_func(x5,y5))


		plt.subplot(236)
		plt.contour(x6,y6,hz_sol,cmap='cool')
		plt.contour(x6,y6,ex_interface_func(x1,y1))
	
	def plot_surface(self):
		pass
