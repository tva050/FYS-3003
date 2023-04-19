import numpy as np
import scipy as sc
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy import integrate

plt.style.use('ggplot')

""" Import data from file """
height, O, N2, O2 = np.loadtxt(r'Data\MSIS.dat', unpack = True ,skiprows=18, usecols= (0, 1, 2, 3))
wavelength, absorption_cross_section_N2, absorption_cross_section_O2, absorption_cross_section_O = np.loadtxt("Data\phot_abs.dat", skiprows=8, unpack=True)
wavelength_fism, irradiance = np.loadtxt("Data\Fism_daily_hr19990216.dat",  unpack = True, skiprows = 51, max_rows=949, delimiter=',', usecols = (1, 2))


""" Constants and conversion """
R_earth = 6371000 # radius of earth (m)
h = 6.62607004e-34 # Planck's constant (J*s)
c = 299792458 # speed of light (m/s) # radius of the earth, [m]


z0 = 0 # [m]
SZA = 75 *(np.pi/180) # solar zenit angle, [rad]
wavelength_fism = wavelength_fism * 1e-9 # -> [m]

height = height * 1e3 # -> [m]
density =  N2*1e6 , O2*1e6, O*1e6 # -> [m^-3]

def _task1_():
    """ 
    Calculate the optical depth as function of wavelength and height
    """
    global unit_optical_depth
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
    def unit_optical_depth(wavelength, height, optical_depth):
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
        plt.plot(wavelength*1e9, unit_optical_depth(wavelength, height, Optical_depth) , color='black', linewidth=1)
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
        plt.plot(density[0], height*1e-3, label='N$_2$')
        plt.plot(density[1], height*1e-3, label='O$_2$', color = "black")
        plt.plot(density[2], height*1e-3, label='O', color = "teal")
        plt.xscale('log')
        plt.xlabel('Density [m$^{-3}$]')
        plt.ylabel('Height [km]')
        plt.title("Thermospheric density with altitude variation")
        plt.legend(loc='upper right', fontsize=10, frameon=True, ncol=3, fancybox=True, shadow=True)
        plt.show()
    return  plot_optical_depth(wavelength, height), plot_absorption_cross_section(), plot_thermospheric_density()

def _task2_(SZA):
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
    
    def unit_optical_depth(wavelength, height, optical_depth):
        def find_nearest(array, value):
            array = np.array(array)
            idx = (np.abs(array - value)).argmin()
            return idx
        tau = []
        for i in range(len(wavelength)):
                index = find_nearest(optical_depth[:,i], 1)
                tau.append(height[index]*1e-3)
        return tau
    # 2: Interpolate the cross section from wavelengths in 'phot_abs.dat' to wavelengths from 'Fism_daily_hr19990216.dat'
    interpolate_N2 = sc.interpolate.interp1d(wavelength, absorption_cross_section_N2, kind='linear')
    interpolate_O2 = sc.interpolate.interp1d(wavelength, absorption_cross_section_O2, kind='linear')
    interpolate_O  = sc.interpolate.interp1d(wavelength, absorption_cross_section_O,  kind='linear')
    # 3: Calculate the optical depth at the interpolated cross section for a specific height and SZAs
    optical_depth_N2 = optical_depth(density[0], interpolate_N2(wavelength_fism))
    optical_depth_O2 = optical_depth(density[1], interpolate_O2(wavelength_fism))
    optical_depth_O  = optical_depth(density[2],  interpolate_O(wavelength_fism))
    interpolated_optical_depth = optical_depth_N2 + optical_depth_O2 + optical_depth_O
    # 4: Calculate the EUV photon flux as function of wavelength and height
    I = irradiance  * np.exp(-interpolated_optical_depth)
    I = (I * (wavelength_fism) / (h * c)) * 1e-4   # W/(nm * m^2) -> Photons/(s * cm^2)
    
    def photon_flux(z):
        I  = np.zeros((len(height),len(wavelength_fism))) 
        for z in range(len(height)):
            I[z,:] = irradiance * np.exp(-interpolated_optical_depth[z,:]) 
            I[z,:] = ((I[z,:] * wavelength_fism) / (h * c)) * 1e-4
        return I 
    
    

    """ plt.pcolormesh(wavelength_fism*1e9, height*1e-3, photon_flux(height))
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Height [km]')
    plt.colorbar(label='Irradiance [$Photons/(s \cdot cm^2$)]')
    plt.title("EUV photon flux, SZA = "+ str(round(SZA*(180/np.pi))) + R"$^\circ$")
    plt.show()
    
    # Photon flux as function of wavelength and height, with unit optical depth
    plt.pcolormesh(wavelength_fism*1e9, height*1e-3, I, cmap='viridis', shading = 'gouraud', norm=colors.LogNorm(vmin=1e8, vmax=1e10))
    plt.plot(wavelength_fism*1e9, unit_optical_depth(wavelength_fism, height ,interpolated_optical_depth), color = 'black', linewidth = 1)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Height (km)')
    plt.colorbar(label = 'Irradiance [$Photons/(s \cdot cm^2$)]')
    plt.title('Photon flux, SZA = ' + str(int(round(np.degrees(SZA), 0))) + ' degrees')
    plt.show()
    
    # Photon flux at different heights
    plt.plot(wavelength_fism*1e9, I[600], label = '600 km')
    plt.plot(wavelength_fism*1e9, I[400], label = '400 km', color = 'teal')
    plt.plot(wavelength_fism*1e9, I[200], label = '200 km', color = 'black')
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('EUV flux (Photons/(s * cm^2))')
    plt.title('Photon flux, SZA = ' + str(int(round(np.degrees(SZA), 0))) + ' degrees')
    plt.show() """
    
    
    np.savetxt("I_SZA75.csv", I[range(0,601)], delimiter=",") #<- Save the photon flux at different heights to a csv file
    
if __name__ == "__main__":
    #_task1_()
    _task2_(SZA)
    