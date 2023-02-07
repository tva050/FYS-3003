import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from matplotlib import colors

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
density = O*1e3 , N2*1e3 , O2*1e3 # -> [m^-3]

def optical_depth(wavelength, height, absorption_cross_section, density):
    tau = np.zeros((len(height),len(wavelength)), dtype=float)
    for i in range(len(height)):
        for j in range(len(wavelength)):
            z0 = height[i]
            integrand = density[height>=z0] * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
            integral = sc.integrate.trapz(integrand, height[height>=z0])
            tau[i][j] = integral * absorption_cross_section[j]
    return tau

sums_up = optical_depth(wavelength, height, absorption_cross_section_col_N2, N2) + optical_depth(wavelength, height, absorption_cross_section_col_O2, O2) + optical_depth(wavelength, height, absorption_cross_section_col_O, O)



plt.pcolormesh(wavelength*1e9, height*1e-3, sums_up , norm=colors.LogNorm(), cmap='viridis')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Height [km]')
plt.title('Optical depth')
plt.colorbar()
plt.show()

""" b) """