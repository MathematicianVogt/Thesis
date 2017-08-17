import numpy as np
import math
import matplotlib.pyplot as plt
import time
class compare_exact_numerical_2D:
	#assumed sol1,sol2 are 2d numpy arrays in i,j format
	def __init__(self,sol1,sol2,xx,yy):
		self.sol1=sol1
		self.sol2=sol2
		self.xx=xx
		self.yy=yy

	def absolute_diff(self):
		dim  = np.shape(self.sol1)
		xdim=dim[0]
		ydim=dim[1]
		print xdim
		print ydim
		res = np.zeros((xdim,ydim))
		for i in range(0,xdim):
			for j in range(0,ydim):
				res[i,j] = math.fabs(self.sol1[i,j]-self.sol2[i,j])
		return res
	def plot(self,cmap='cool'):
		plt.figure(1)
		plt.subplot(121)
		plt.pcolor(self.xx,self.yy,self.sol1)
		plt.xlabel("x")
		plt.ylabel("y")
		plt.title("Numerical Solution")

		plt.subplot(122)
		plt.pcolor(self.xx,self.yy,self.sol2)
		plt.xlabel("x")
		plt.ylabel("y")
		plt.title("Exact Solution")

		plt.show()







