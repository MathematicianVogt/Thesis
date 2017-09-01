import numpy as np
import math
import matplotlib.pyplot as plt
import time
class compare_exact_numerical_2D:
	#assumed sol1 are 2d numpy arrays in i,j format. Exact sol is a lambda function
	def __init__(self,sol1,exact_sol,xx,yy,x,y):
		self.sol1=sol1
		self.exact=exact_sol
		self.exact_sol = self.construct_exact(self.exact)
		self.xx=xx
		self.yy=yy
		self.x_list=x
		self.y_list = y

	
	def construct_exact(self,function_ev):
		dim = np.shape(self.sol1)
		xdim=dim[0]
		ydim=dim[1]
		print xdim
		print ydim
		res = np.zeros((xdim,ydim))
		for i in range(0,xdim):
			for j in range(0,ydim):
				res[i,j] = function_ev( self.x_list[i],self.y_list[j])
		return res



	def absolute_diff(self):
		dim  = np.shape(self.sol1)
		xdim=dim[0]
		ydim=dim[1]
		print xdim
		print ydim
		res = np.zeros((xdim,ydim))
		for i in range(0,xdim):
			for j in range(0,ydim):
				res[i,j] = math.fabs(self.sol1[i,j]-self.exact_sol[i,j])
		return res
	def plot(self,cmap='cool'):
		plt.figure(1)
		plt.subplot(131)
		plt.pcolor(self.xx,self.yy,self.sol1)
		plt.xlabel("x")
		plt.ylabel("y")
		plt.title("Numerical Solution")

		plt.subplot(132)
		plt.pcolor(self.xx,self.yy,self.exact_sol)
		plt.xlabel("x")
		plt.ylabel("y")
		plt.title("Exact Solution")

		plt.subplot(133)
		plt.pcolor(self.xx,self.yy,self.absolute_diff())
		plt.xlabel("x")
		plt.ylabel("y")
		plt.title("Absolute Error")

		plt.show()







