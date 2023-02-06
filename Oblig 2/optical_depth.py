import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

plt.style.use('ggplot')



""" Import data from file """
height, O, N2, O2, mass_density, temperature_neutral = np.loadtxt('MSIS.dat', skiprows=18, unpack=True)
wavelength, absorption_cross_section_col_N2, absorption_cross_section_col_O2, absorption_cross_section_col_O = np.loadtxt("phot_abs.dat", skiprows=8, unpack=True)
sigma_ionization = sc.io.loadmat('sigma_ionization.mat')

""" Constants and conversion """
R_earth = 6.371e6 # radius of the earth, [m]

z0 = 0 # [m]
SZA = 0 # solar zenit angles

height = height * 1e3 # -> [m]
number_density = (O + N2 + O2)*1e-6 # -> [m^-3]
#wavelength = wavelength * 1e-1 # -> [nm]

# Expression for the optical depth when the solar zenith angle is smaller than 90 degrees

integrand = mass_density[height>=z0] * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
integral = sc.integrate.trapz(integrand, height[height>=z0])


multiplied = integral * absorption_cross_section_col_N2

def optical_depth(wavelength, height, absorption_cross_section):
    tau = np.zeros((len(height), len(wavelength)), dtype=float)
    for z0 in range(len(height)):
        for j in range(len(wavelength)):
            integrand = number_density[height>=z0] * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
            integral = sc.integrate.trapz(integrand, height[height>=z0])
            tau[z0][j] = integral * absorption_cross_section[j]
    return tau


sums_up = optical_depth(wavelength, height, absorption_cross_section_col_N2) + optical_depth(wavelength, height, absorption_cross_section_col_O2) + optical_depth(wavelength, height, absorption_cross_section_col_O)

plt.pcolormesh(wavelength*1e9, height*1e-3, sums_up)
plt.colorbar()
plt.title('Optical depth')
plt.ylabel('Height [km]')
plt.xlabel('Wavelength [nm]')
plt.show()
 
 
""" # Plot the absorption cross section of N2 with varying wavelength from "phot_abs.dat"
plt.plot(wavelength, multiplied*1e-3)
plt.title('Absorption cross section of N2 with varying wavelength')
plt.ylabel(R"$N_2$ optical depth [$km^{-1}$]")
plt.xlabel('Wavelength [Å]')
plt.xscale('log')
plt.yscale('log')
plt.show()
 """


""" b) """

