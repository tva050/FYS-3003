import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from matplotlib import colors

plt.style.use('ggplot')

""" Import data from file """
height, O, N2, O2, mass_density, temperature_neutral = np.loadtxt('Data\MSIS.dat', skiprows=18, unpack=True)
wavelength, absorption_cross_section_col_N2, absorption_cross_section_col_O2, absorption_cross_section_col_O = np.loadtxt("Data\phot_abs.dat", skiprows=8, unpack=True)
time, wavelength_fism, irradiance, uncertainty  = np.loadtxt("Data\Fism_daily_hr19990216.dat", delimiter="," ,skiprows= 1,unpack=True)
sigma_ionization = sc.io.loadmat('Data\sigma_ionization.mat')

""" Constants and conversion """
R_earth = 6.371e6 # radius of the earth, [m]

z0 = 0 # [m]
SZA = 0 *(np.pi/180) # solar zenit angle, [rad]

height = height * 1e3 # -> [m]
density = O*1e6 , N2*1e6 , O2*1e6 # -> [m^-3]
""" a) """
def optical_depth(wavelength, height, absorption_cross_section, density):
    tau = np.zeros((len(height),len(wavelength)), dtype=float)
    for i in range(len(height)):
        for j in range(len(wavelength)):
            z0 = height[i]
            integrand = density[height>=z0] * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
            integral = sc.integrate.trapz(integrand, height[height>=z0])
            tau[i][j] = integral * absorption_cross_section[j]
    return tau

#sums_up = optical_depth(wavelength, height, absorption_cross_section_col_N2, density[0]) + optical_depth(wavelength, height, absorption_cross_section_col_O2, density[1]) + optical_depth(wavelength, height, absorption_cross_section_col_O, density[2])
def get_optical_depth(wavelength, height, absorption_cs1, absorption_cs2, absorption_cs3):
    get_od = optical_depth(wavelength, height, absorption_cs1, density[0]) + optical_depth(wavelength, height, absorption_cs2, density[1]) + optical_depth(wavelength, height, absorption_cs3, density[2])
    return get_od

# plot a line at unit optical depth (i.e. where tau is closest to one) for each wavelength
def tau_close_to_one(wavelength, height, optical_depth):
    def find_nearest(array, value):
        array = np.array(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    tau = []
    for i in range(len(wavelength)):
            index = find_nearest(optical_depth[:,i], 1)
            tau.append(height[index]*1e-3)
    return tau

def plot_optical_depth(wavelength, height):
    plt.pcolormesh(wavelength*1e9, height*1e-3, get_optical_depth(wavelength, height ,absorption_cross_section_col_N2, absorption_cross_section_col_O2, absorption_cross_section_col_O) , norm=colors.LogNorm())
    plt.plot(wavelength*1e9, tau_close_to_one(wavelength, height, get_optical_depth(wavelength, height ,absorption_cross_section_col_N2, absorption_cross_section_col_O2, absorption_cross_section_col_O)), color='black', linewidth=1)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Height [km]')
    plt.colorbar(label='Optical depth')
    plt.xlim(4.9,100)
    plt.title("Optical depth, SZA = "+ str(round(SZA*(180/np.pi))) + R"$^\circ$")
    plt.show()

def plot_absorption_cross_section():
    plt.plot(wavelength*1e9, absorption_cross_section_col_N2, label='N$_2$')
    plt.plot(wavelength*1e9, absorption_cross_section_col_O2, label='O$_2$', color = "black")
    plt.plot(wavelength*1e9, absorption_cross_section_col_O, label='O', color = "teal")
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Absorption cross section [m$^2$]')
    plt.title("Absorption cross section")
    plt.legend(loc='upper left', fontsize=10, frameon=True, ncol=3, fancybox=True, shadow=True)
    plt.show()
    
def plot_thermospheric_density():
    plt.plot(density[1], height*1e-3, label='N$_2$')
    plt.plot(density[2], height*1e-3, label='O$_2$', color = "black")
    plt.plot(density[0], height*1e-3, label='O', color = "teal")
    plt.xscale('log')
    plt.xlabel('Density [m$^{-3}$]')
    plt.ylabel('Height [km]')
    plt.title("Thermospheric density with altitude variation")
    plt.legend(loc='upper right', fontsize=10, frameon=True, ncol=3, fancybox=True, shadow=True)
    plt.show()


""" 
b)
- Calculate the the EUV-photon flux as function of wavelength and height 
(The photon flux is given from the thermo-sphere to infinity, i.e. from 80 km to infinity)
"""

# 1: extracting the data from fism file [DONE]

# 2: Interpolate the cross section from wavelengths in 'phot_abs.dat' to wavelengths from 'Fism_daily_hr19990216.dat'

def interpolation(wavelength, absorption_cross_section):
    ip = sc.interpolate.interp1d(wavelength, absorption_cross_section, kind='linear') 
    return ip


# 3: Calculate the optical depth at the interpolated cross section for a specific height and SZA

# 4: Calculate the EUV photon flux as function of wavelength and height

#def EUV_photon_flux():

# 5: 5.	You can use the relationship between the energy of a photon and its frequency or wavelength to convert the flux from W/(nm * m^2) to Photons/(s * cm^2).