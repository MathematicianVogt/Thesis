import numpy as np
import math

'''
get_theta:
Inputs: x_prime (x'(t)), y_prime(y'(t))
Output: theta array
Description: Find the theta value at the interface, for each point on the interface. x_prime, y_prime can be scalars of arrays, i.e can find one theta, or 
all theta for interface  


Can get theta at every point on interface, or point by point, based on if x_prime, y_prime are scalars, or the entire set x(t), y(t)
'''
def get_theta(x_prime,y_prime):
	cos_theta = np.divide(y_prime,np.power(np.power(x_prime,2)+np.power(y_prime,2),.5) +math.pow(10.0,-14))
	#sin_theta = np.divide(-x_prime ,np.power(np.power(x_prime,2)+np.power(y_prime,2),.5) +math.pow(10.0,-14))
	theta = np.arccos(cos_theta)
	return theta


'''
get_xi:
Inputs:
	x - x cor near interface
	y - y cor near interface
	X - X cor on interface
	Y- Y cor on interface
	theta - theta value on interface
Output:
	xi - the local tranformation xi for the local point (x,y), near the interface (X,Y)
'''
def get_xi(x,y,X,Y,theta):
	xi = np.add(np.prod(x-X,np.cos(theta)),np.prod(y-Y,np.sin(theta)))
	return xi

'''
get_eta:
Inputs:
	x - x cor near interface
	y - y cor near interface
	X - X cor on interface
	Y- Y cor on interface
	theta - theta value on interface
Output:
	eta - the local tranformation eta for the local point (x,y), near the interface (X,Y)
'''	
def get_eta(x,y,X,Y,theta):
	eta = np.add(np.prod(-(x-X),np.sin(theta)),np.prod(y-Y,np.cos(theta)))
	return eta







