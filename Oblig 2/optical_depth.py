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

""" a) """

z0 = 0 # m
SZA = 0 # solar zenit angles, degrees

height = height * 1e3 # conversion to m
#mass_density = mass_density * 1e-3  # conversion to kg/m^3
number_density = (O + N2 + O2)*1e-6

# Expression for the optical depth when the solar zenith angle is smaller than 90 degrees
integrand = mass_density[height>=z0] * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
integral = sc.integrate.trapz(integrand, height[height>=z0])


multiplied = integral * absorption_cross_section_col_N2

wavelength = wavelength * 1e10 # conversion to angstrom (Å)

def optical_depth(z0, SZA):
    
    
# Plot the absorption cross section of N2 with varying wavelength from "phot_abs.dat"
plt.plot(wavelength, absorption_cross_section_col_N2)
plt.title('Absorption cross section of N2 with varying wavelength')
plt.ylabel(R"$N_2$ cross section [m$^2$]")
plt.xlabel('Wavelength [Å]')
#plt.xscale('log')
plt.yscale('log')
plt.show()

plt.plot(wavelength, multiplied)
plt.ylabel(R"$\tau$")
plt.yscale('log')
plt.show()



""" b) """
