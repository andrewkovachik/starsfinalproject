
import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt


sigma = 5.6703*10**-8

# radius_surface = 6900000000
# temperature_surface = 56000 
# energygen_centre = 0.001
# x = np.linspace(0, 100000, 10000)
# y = [0.0]*10000
# for i in range(len(x)):
# 	y[i] = 	(4/3)*np.pi*(radius_surface**3)*x[i]*energygen_centre - (4.0*np.pi*sigma*(radius_surface**2)*(temperature_surface**4))/(np.sqrt(4.0*np.pi*sigma*(radius_surface**2)*(temperature_surface**4)*(4/3)*np.pi*(radius_surface**3)*i*energygen_centre))
# plt.plot(x,y)
# plt.show()


# Surface Temperature = 
# Central Density = 74237.74616514979
# Radius = 562477594.72664583
# Mass = 0.70906952065099615
# Luminosity = 4.6751081103682353e+25



def Lum_error(density_centre):

	radius_surface = 562477594.72664583
	temperature_surface = 3794.8714590341042
	energygen_centre = 0.001

	L_star = (4/3)*np.pi*(radius_surface**3)*density_centre*energygen_centre
	return L_star - (4.0*np.pi*sigma*(radius_surface**2)*(temperature_surface**4))/(np.sqrt(4.0*np.pi*sigma*(radius_surface**2)*(temperature_surface**4)*L_star))

def Bisection(function, initial_guess):

	print(optimize.bisect(function, initial_guess, initial_guess + (initial_guess/10)))


#print(optimize.bisect(Lum_error(), 10.0**12, 6500 + (10.0**12/2)))
	

# Lum_error(69000000,5700,40000,10**12)
Bisection(Lum_error, 74237.74616514979)