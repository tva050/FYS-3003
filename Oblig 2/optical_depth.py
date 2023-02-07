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
density = O*1e6 , N2*1e6 , O2*1e6 # -> [m^-3]
""" a) """

def optical_depth(wavelength, height, absorption_cross_section, density):
    tau = np.zeros((len(height),len(wavelength)), dtype=float)
    for i in range(len(height)):
        for j in range(len(wavelength)):
            z0 = height[i]
            integrand = density[height>=z0]*(1e6) * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
            integral = sc.integrate.trapz(integrand, height[height>=z0])
            tau[i][j] = integral * absorption_cross_section[j]
    return tau

sums_up = optical_depth(wavelength, height, absorption_cross_section_col_N2, N2) + optical_depth(wavelength, height, absorption_cross_section_col_O2, O2) + optical_depth(wavelength, height, absorption_cross_section_col_O, O)

# plot a line at unit optical depth (i.e. where tau is closest to one) for each wavelength
def tau_close_to_one(wavelength, height, sums_up):
    def find_nearest(array, value):
        array = np.array(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    tau = []
    for i in range(len(wavelength)):
            index = find_nearest(sums_up[:,i], 1)
            tau.append(height[index]*1e-3)
    return tau

plt.pcolormesh(wavelength*1e9, height*1e-3, sums_up , norm=colors.LogNorm())
plt.plot(wavelength*1e9, tau_close_to_one(wavelength, height, sums_up), color='black')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Height [km]')
plt.colorbar(label='Optical depth')
plt.xlim(4.9,100)
plt.show()

""" b) """