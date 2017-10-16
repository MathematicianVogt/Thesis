import pylab as plt
import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib.animation as manimation
import numpy as np
import matplotlib
import math
class maxwell_plot:
	def __init__(self,ex,ey,ez,hx,hy,hz,Tmax,nt):
		self.ex=ex
		self.ey=ey
		self.ez=ez
		self.hx=hx
		self.hy=hy
		self.hz=hz
		self.t_list = np.linspace(0,Tmax,nt)
		self.t_list_shift = [x + (self.t_list[1]-self.t_list[0])/2.0 for x in self.t_list]
		self.t_list_shift = self.t_list_shift[:-1]
		


	def plot_contour(self,filename):
		default_cmap='hot'
		print "Starting to Plot Contour Plot All....."
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

		(x6,y6,hz_sol) = self.hz.get_sol()
		hz_interface_func = self.hz.get_phi()
		

		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=filename, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=60, metadata=metadata)
		fig=plt.figure()
		with writer.saving(fig, filename+ ".mp4", 100):
			


			for i in range(0,len(self.t_list)-1):
				print "ploting .... t = " + str(self.t_list[i])
				fulltime=self.t_list[i+1]
				halftime=self.t_list_shift[i]
				
				plt.subplot(231)
				plt.pcolor(x1,y1,ex_sol[i],cmap=default_cmap)
				plt.colorbar()

				#plt.contour(x1,y1,ex_interface_func(x1,y1),[1])
				plt.title(r"$E_x(x,y,t) : t=$" + str(fulltime))

				plt.subplot(232)
				plt.pcolor(x2,y2,ey_sol[i],cmap=default_cmap)
				plt.colorbar()
				#plt.contour(x2,y2,ey_interface_func(x2,y2),[1])
				plt.title(r"$E_y(x,y,t) : t=$" + str(fulltime))

				plt.subplot(233)
				plt.pcolor(x3,y3,ez_sol[i],cmap=default_cmap)
				plt.colorbar()
				#plt.contour(x3,y3,ez_interface_func(x3,y3),[1])
				plt.title(r"$E_z(x,y,t) : t=$" + str(fulltime))


				plt.subplot(234)
				plt.pcolor(x4,y4,hx_sol[i],cmap=default_cmap)
				plt.colorbar()
				#plt.contour(x4,y4,hx_interface_func(x4,y4),[1])
				plt.title(r"$H_x(x,y,t) : t=$" + str(halftime))


				plt.subplot(235)
				plt.pcolor(x5,y5,hy_sol[i],cmap=default_cmap)
				plt.colorbar()

				#plt.contour(x5,y5,hy_interface_func(x5,y5),[1])
				plt.title(r"$H_y(x,y,t) : t=$" + str(halftime))


				plt.subplot(236)
				plt.pcolor(x6,y6,hz_sol[i],cmap=default_cmap)
				plt.colorbar()
				#plt.contour(x6,y6,ex_interface_func(x1,y1),[1])
				plt.title(r"$H_z(x,y,t) : t=$" + str(halftime))
				writer.grab_frame()
				plt.clf()
	
	def plot_contour_TE(self,filename):
		default_cmap='hot'
		print "Starting to Plot Contour Plot TE Solutions....."
		(x1,y1,ex_sol) = self.ex.get_sol()
		ex_interface_func = self.ex.get_phi()

		(x2,y2,ey_sol) = self.ey.get_sol()
		ey_interface_func = self.ey.get_phi()

		(x5,y5,hz_sol) = self.hz.get_sol()
		hz_interface_func = self.hz.get_phi()
		

		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=filename, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=60, metadata=metadata)
		

		fig=plt.figure()
		z_min=-9.0
		z_max=13.0

		with writer.saving(fig, filename+ ".mp4", 100):
			
			plt.subplot(131)
			plt.pcolor(x5,y5,hz_sol[0],cmap=default_cmap,vmin=z_min,vmax=z_max)
			#plt.contour(x5,y5,hz_interface_func(x5,y5),[1])
			plt.colorbar()
			
			plt.title(r"$H_z(x,y,t) : t=$0")

			plt.subplot(132)
			plt.pcolor(x1,y1,ex_sol[0],cmap=default_cmap,vmin=z_min,vmax=z_max)
			#plt.contour(x1,y1,ey_interface_func(x1,y1),[1])
			plt.colorbar()
			
			plt.title(r"$E_x(x,y,t) : t=$0")

			plt.subplot(133)
			plt.pcolor(x2,y2,ey_sol[0],cmap=default_cmap,vmin=z_min,vmax=z_max)
			plt.colorbar()
			
			#plt.contour(x3,y3,ey_interface_func(x3,y3),[1])
			plt.title(r"$E_y(x,y,t) : t=$0" )
			writer.grab_frame()
			plt.clf()


			for i in range(1,len(self.t_list)):
				print "ploting .... t = " + str(self.t_list[i])
				
				fulltime=self.t_list[i]
				halftime=self.t_list_shift[i-1]
				
				plt.subplot(131)
				plt.pcolor(x5,y5,hz_sol[i],cmap=default_cmap,vmin=z_min,vmax=z_max)
				#plt.contour(x5,y5,hz_interface_func(x5,y5),[1])
				plt.colorbar()
				
				plt.title(r"$H_z(x,y,t) : t=$" + str(halftime))

				plt.subplot(132)
				plt.pcolor(x1,y1,ex_sol[i],cmap=default_cmap,vmin=z_min,vmax=z_max)
				#plt.contour(x1,y1,ey_interface_func(x1,y1),[1])
				plt.colorbar()
				
				plt.title(r"$E_x(x,y,t) : t=$" + str(fulltime))

				plt.subplot(133)
				plt.pcolor(x2,y2,ey_sol[i],cmap=default_cmap,vmin=z_min,vmax=z_max)
				plt.colorbar()
				
				#plt.contour(x3,y3,ey_interface_func(x3,y3),[1])
				plt.title(r"$E_y(x,y,t) : t=$" + str(fulltime))

				writer.grab_frame()
				plt.clf()

	def plot_contour_TM(self,filename):
		default_cmap='hot'
		print "Starting to Plot Contour Plot TE Solutions....."
		(x3,y3,ez_sol) = self.ez.get_sol()
		ez_interface_func = self.ez.get_phi()

		

		(x4,y4,hx_sol) = self.hx.get_sol()
		hx_interface_func = self.hx.get_phi()

		(x5,y5,hy_sol) = self.hy.get_sol()
		hy_interface_func = self.hy.get_phi()


		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=filename, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=60, metadata=metadata)
		fig=plt.figure()
		with writer.saving(fig, filename+ ".mp4", 100):
			for i in range(0,len(self.t_list)):
				print "ploting .... t = " + str(self.t_list[i])
				
				fulltime=self.t_list[i]
				halftime=self.t_list_shift[i]
				
				plt.subplot(231)
				plt.pcolor(x4,y4,hx_sol[i],cmap=default_cmap)
				#plt.contour(x1,y1,ex_interface_func(x1,y1),[1])
				plt.title(r"$H_x(x,y,t) : t=$" + str(fulltime))
				plt.colorbar()

				plt.subplot(232)
				plt.pcolor(x5,y5,hy_sol[i],cmap=default_cmap)
				#plt.contour(x2,y2,ey_interface_func(x2,y2),[1])
				plt.title(r"$H_y(x,y,t) : t=$" + str(fulltime))
				plt.colorbar()

				plt.subplot(233)
				plt.pcolor(x3,y3,ez_sol[i],cmap=default_cmap)
				#plt.contour(x3,y3,ez_interface_func(x3,y3),[1])
				plt.title(r"$E_z(x,y,t) : t=$" + str(fulltime))
				plt.colorbar()

				writer.grab_frame()
				plt.clf()
	


	def plot_surface(self):
		default_cmap='hot'
		print "Starting to Plot All Surfaces....."
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

		(x6,y6,hz_sol) = self.hz.get_sol()
		hz_interface_func = self.hz.get_phi()
		

		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=filename, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=60, metadata=metadata)
		fig=plt.figure()
		with writer.saving(fig, filename+ ".mp4", 100):
			for i in range(0,len(self.t_list)):
				print "ploting .... t = " + str(self.t_list[i])
				fulltime=self.t_list[i+1]
				halftime=self.t_list_shift[i]
				
				plt.subplot(231)
				plt.plot_surface(x1,y1,ex_sol[i],cmap=default_cmap)
				

				#plt.contour(x1,y1,ex_interface_func(x1,y1),[1])
				plt.title(r"$E_x(x,y,t) : t=$" + str(fulltime))

				plt.subplot(232)
				plt.plot_surface(x2,y2,ey_sol[i],cmap=default_cmap)
				
				#plt.contour(x2,y2,ey_interface_func(x2,y2),[1])
				plt.title(r"$E_y(x,y,t) : t=$" + str(fulltime))

				plt.subplot(233)
				plt.plot_surface(x3,y3,ez_sol[i],cmap=default_cmap)
				
				#plt.contour(x3,y3,ez_interface_func(x3,y3),[1])
				plt.title(r"$E_z(x,y,t) : t=$" + str(fulltime))


				plt.subplot(234)
				plt.plot_surface(x4,y4,hx_sol[i],cmap=default_cmap)
				
				#plt.contour(x4,y4,hx_interface_func(x4,y4),[1])
				plt.title(r"$H_x(x,y,t) : t=$" + str(halftime))


				plt.subplot(235)
				plt.plot_surface(x5,y5,hy_sol[i],cmap=default_cmap)
				

				#plt.contour(x5,y5,hy_interface_func(x5,y5),[1])
				plt.title(r"$H_y(x,y,t) : t=$" + str(halftime))


				plt.subplot(236)
				plt.plot_surface(x6,y6,hz_sol[i],cmap=default_cmap)
				
				#plt.contour(x6,y6,ex_interface_func(x1,y1),[1])
				plt.title(r"$H_z(x,y,t) : t=$" + str(halftime))
				writer.grab_frame()
				plt.clf()


	def plot_surface_TE(self,filename):
			default_cmap='hot'
			print "Starting to Plot Surface Plot TE Solutions....."
			(x1,y1,ex_sol) = self.ex.get_sol()
			ex_interface_func = self.ex.get_phi()

			(x2,y2,ey_sol) = self.ey.get_sol()
			ey_interface_func = self.ey.get_phi()

			(x5,y5,hz_sol) = self.hz.get_sol()
			hz_interface_func = self.hz.get_phi()
			

			FFMpegWriter = manimation.writers['ffmpeg']
			metadata = dict(title=title, artist='Matplotlib',
			                comment='Movie support!')
			writer = FFMpegWriter(fps=60, metadata=metadata)
			fig=plt.figure()
			with writer.saving(fig, filename+ ".mp4", 100):
				for i in range(0,len(self.t_list)):
					
					fulltime=self.t_list[i]
					halftime=self.t_list_shift[i]
					
					plt.subplot(131)
					plt.plot_surface(x5,y5,hz_sol[i],cmap=default_cmap)
					#plt.contour(x5,y5,hz_interface_func(x5,y5),[1])
					plt.title(r"$H_z(x,y,t) : t=$" + str(halftime))

					plt.subplot(132)
					plt.plot_surface(x1,y1,ex_sol[i],cmap=default_cmap)
					#plt.contour(x1,y1,ey_interface_func(x1,y1),[1])
					plt.title(r"$E_x(x,y,t) : t=$" + str(fulltime))

					plt.subplot(133)
					plt.plot_surface(x3,y3,ey_sol[i],cmap=default_cmap)
					#plt.contour(x3,y3,ey_interface_func(x3,y3),[1])
					plt.title(r"$E_y(x,y,t) : t=$" + str(fulltime))

					writer.grab_frame()
					plt.clf()

	def plot_surface_TM(self,filename):
		default_cmap='hot'
		print "Starting to Plot Surface Plot TM Solutions....."
		(x3,y3,ez_sol) = self.ez.get_sol()
		ez_interface_func = self.ez.get_phi()

		

		(x4,y4,hx_sol) = self.hx.get_sol()
		hx_interface_func = self.hx.get_phi()

		(x5,y5,hy_sol) = self.hy.get_sol()
		hy_interface_func = self.hy.get_phi()


		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=title, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=60, metadata=metadata)
		fig=plt.figure()
		with writer.saving(fig, filename+ ".mp4", 100):
			for i in range(0,len(self.t_list)):
				
				fulltime=self.t_list[i]
				halftime=self.t_list_shift[i]
				
				plt.subplot(231)
				plt.plot_surface(x4,y4,hx_sol[i],cmap=default_cmap)
				#plt.contour(x1,y1,ex_interface_func(x1,y1),[1])
				plt.title(r"$H_x(x,y,t) : t=$" + str(fulltime))

				plt.subplot(232)
				plt.plot_surface(x5,y5,hy_sol[i],cmap=default_cmap)
				#plt.contour(x2,y2,ey_interface_func(x2,y2),[1])
				plt.title(r"$H_y(x,y,t) : t=$" + str(fulltime))

				plt.subplot(233)
				plt.plot_surface(x3,y3,ez_sol[i],cmap=default_cmap)
				#plt.contour(x3,y3,ez_interface_func(x3,y3),[1])
				plt.title(r"$E_z(x,y,t) : t=$" + str(fulltime))

				writer.grab_frame()
				plt.clf()
	
