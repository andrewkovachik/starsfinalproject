
import scipy.optimize as optimize
import numpy as np

def Lum_error(radius_surface,temperature_surface,energygen_centre, density_centre):

	L_star = (4/3)*np.pi*(radius_surface**3)*density_centre*energygen_centre
	return L_star - (4*np.pi*sigma*(radius_surface**2)*(temperature_surface**4))/(np.sqrt(4*np.pi*sigma*(radius_surface**2)*(temperature_surface**4)*L_star))

def Bisect(function, initial_guess)

	print(optimize.bisect(function, initial_guess, initial_guess + initial_guess/2))


		