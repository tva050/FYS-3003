import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

plt.style.use('ggplot')

""" Constants """

R_earth = 6371 # km
kji = np.arange(0, 90, 15) # degrees

""" Import data from file """
# Data from file MSIS.dat 
height, O, N2, O2, mass_density, temperature_neutral = np.loadtxt('MSIS.dat', skiprows=18, unpack=True)
# Data from file phot_abs.dat
wavelength, absorption_cross_section_col2, absorption_cross_section_col3, absorption_cross_section_col4 = np.loadtxt("phot_abs.dat", skiprows=8, unpack=True)
# Data from file sigma_ionization.mat
sigma_ionization = sp.io.loadmat('sigma_ionization.mat')


""" 

a) Make functions that calculate the optical depth as a function of altitude and wavelength,
for vertical incidence, and for variable zenith-angle of the incident light. 

"""

# phi are the wavelength-dependent absorption cross sections of species j
# n_j(0) are the heights profile of the concentration of species j
# kji is the zenith angle of the incident light
# z is the altitude 



    


