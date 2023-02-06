import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

plt.style.use('ggplot')

""" Constants """

# radius of the earth in m
R_earth = 6.371e6 # m

""" Import data from file """
# Data from file MSIS.dat 
height, O, N2, O2, mass_density, temperature_neutral = np.loadtxt('MSIS.dat', skiprows=18, unpack=True)
# Data from file phot_abs.dat
wavelength, absorption_cross_section_col_N2, absorption_cross_section_col_O2, absorption_cross_section_col_O = np.loadtxt("phot_abs.dat", skiprows=8, unpack=True)
# Data from file sigma_ionization.mat
sigma_ionization = sc.io.loadmat('sigma_ionization.mat')


""" 
a) Make functions that calculate the optical depth as a function of altitude and wavelength,
for vertical incidence, and for variable zenith-angle of the incident light. 

- phi: are the wavelength-dependent absorption cross sections of species j
- n_j(0): are the heights profile of the concentration of species j
- kji: is the zenith angle of the incident light
- z: is the altitude 
"""

z0 = 75_000 # m
SZA = 60 # solar zenit angles, degrees

# heigth conversion to m
height = height * 1e3
# mass density conversion to kg/m^3
#mass_density = mass_density * 1e-3
number_density = (O + N2 + O2)*1e6

integrand = number_density[height>=z0] * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
integral = sc.integrate.trapz(integrand, height[height>=z0])

# Now you just multiply this integral with the absorption cross section at each wavelength and you'll get the optical depth at a range 
# of wavelengths for a particular SZA and altitude.

multiplied = integral * absorption_cross_section_col_N2

# Wavelegth conversion to angstroms
wavelength = wavelength * 1e10

plt.plot(wavelength, multiplied)
plt.title('Absorption cross section of N2 with varying wavelength')
plt.ylabel(R"$N_2$ cross section [m$^2$]")
plt.xlabel('Wavelength [Å]')
#plt.xscale('log')
plt.yscale('log')
plt.show()

""" 
b) 
"""

