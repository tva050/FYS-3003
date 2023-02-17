import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from matplotlib import colors

plt.style.use('ggplot')

""" Import data from file """
height, O, N2, O2, mass_density, temperature_neutral = np.loadtxt('Data\MSIS.dat', skiprows=18, unpack=True)
wavelength, absorption_cross_section_N2, absorption_cross_section_O2, absorption_cross_section_O = np.loadtxt("Data\phot_abs.dat", skiprows=8, unpack=True)
time, wavelength_fism, irradiance, uncertainty  = np.loadtxt("Data\Fism_daily_hr19990216.dat", delimiter="," ,skiprows= 1,unpack=True)
sigma_ionization = sc.io.loadmat('Data\sigma_ionization.mat')

""" Constants and conversion """
R_earth = 6.371e6 # radius of the earth, [m]
h = 6.626e-34 # Planck constant, [J*s]
c = 299_792_458 # speed of light, [m/s]

z0 = 0 # [m]
SZA = 0 *(np.pi/180) # solar zenit angle, [rad]

height = height * 1e3 # -> [m]
density =  N2*1e6 , O2*1e6, O*1e6 # -> [m^-3]
""" a) """

def _task1_():
    def optical_depth(wavelength, height, absorption_cross_section, density):
        tau = np.zeros((len(height),len(wavelength)), dtype=float)
        for i in range(len(height)):
            for j in range(len(wavelength)):
                z0 = height[i]
                integrand = density[height>=z0] * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
                integral = sc.integrate.trapz(integrand, height[height>=z0])
                tau[i][j] = integral * absorption_cross_section[j]
        return tau

    Optical_depth = optical_depth(wavelength, height, absorption_cross_section_N2, density[0]) + optical_depth(wavelength, height, absorption_cross_section_O2, density[1]) + optical_depth(wavelength, height, absorption_cross_section_O, density[2])

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
        plt.pcolormesh(wavelength*1e9, height*1e-3, Optical_depth , norm=colors.LogNorm())
        plt.plot(wavelength*1e9, tau_close_to_one(wavelength, height, Optical_depth) , color='black', linewidth=1)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Height [km]')
        plt.colorbar(label='Optical depth')
        plt.xlim(4.9,100)
        plt.title("Optical depth, SZA = "+ str(round(SZA*(180/np.pi))) + R"$^\circ$")
        plt.show()

    def plot_absorption_cross_section():
        plt.plot(wavelength*1e9, absorption_cross_section_N2, label='N$_2$')
        plt.plot(wavelength*1e9, absorption_cross_section_O2, label='O$_2$', color = "black")
        plt.plot(wavelength*1e9, absorption_cross_section_O, label='O', color = "teal")
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
    return  plot_optical_depth(wavelength, height), plot_absorption_cross_section(), plot_thermospheric_density()

def _task2_():
    """ 
    b)
    - Calculate the the EUV-photon flux as function of wavelength and height 
    (The photon flux is given from the thermo-sphere to infinity, i.e. from 80 km to infinity)
    """
    def optical_depth(density, absorption_cross_section):
        tau = np.zeros((len(height),len(wavelength_fism)), dtype=float)
        for i in range(len(height)):
            for j in range(len(wavelength_fism)):
                z0 = height[i]
                integrand = density[height>=z0] * (1 - ( (R_earth + z0)/(R_earth + height[height>=z0]) ) ** 2 * (np.sin(SZA)) ** 2 ) ** (-0.5)
                integral = sc.integrate.trapz(integrand, height[height>=z0])
                tau[i][j] = integral * absorption_cross_section[j]
        return tau
    # 2: Interpolate the cross section from wavelengths in 'phot_abs.dat' to wavelengths from 'Fism_daily_hr19990216.dat'
    interpolate_N2 = sc.interpolate.interp1d(wavelength, absorption_cross_section_N2, kind='linear', fill_value="extrapolate")
    interpolate_O2 = sc.interpolate.interp1d(wavelength, absorption_cross_section_O2, kind='linear', fill_value="extrapolate")
    interpolate_O  = sc.interpolate.interp1d(wavelength, absorption_cross_section_O,  kind='linear', fill_value="extrapolate")
    # 3: Calculate the optical depth at the interpolated cross section for a specific height and SZAs
    optical_depth_N2 = optical_depth(density[0], interpolate_N2(wavelength_fism))
    optical_depth_O2 = optical_depth(density[1], interpolate_O2(wavelength_fism))
    optical_depth_O  = optical_depth(density[2],  interpolate_O(wavelength_fism))
    #print(optical_depth_N2)
    interpolated_optical_depth = optical_depth_N2 + optical_depth_O2 + optical_depth_O
    #print(interpolated_optical_depth)
    # 4: Calculate the EUV photon flux as function of wavelength and height
    I = irradiance * np.exp(-interpolated_optical_depth)
    #print(I)
    #I = ((I * wavelength_fism) / (h * c)) * 1e-4   # W/(nm * m^2) -> Photons/(s * cm^2)
    #print(I)
    
    """ def photon_flux(z):
        I  = np.zeros((len(height),len(wavelength_fism))) 
        for z in range(len(height)):
            I[z,:] = irradiance * np.exp(-interpolated_optical_depth[z,:]) 
            I[z,:] = ((I[z,:] * wavelength_fism) / (h * c)) * 1e-4
        return I """

    
    plt.pcolormesh(wavelength_fism*1e9, height*1e-3, I, norm=colors.LogNorm())
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Height [km]')
    plt.colorbar(label='EUV photon flux [Photons/(s * cm$^2$)]')
    plt.title("EUV photon flux, SZA = "+ str(round(SZA*(180/np.pi))) + R"$^\circ$")
    plt.show()

if __name__ == "__main__":
    #_task1_()
    _task2_()