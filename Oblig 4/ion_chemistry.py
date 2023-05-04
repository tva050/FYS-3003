import numpy as np
import matplotlib.pyplot as plt

plt.style.use('ggplot')

height, O, N2, O2, mass_density, nutral_temp = np.loadtxt(r'Data\MSIS.dat', unpack = True ,skiprows=118, usecols= (0, 1, 2, 3, 4, 5))
wavelength, absorption_cross_section_N2, absorption_cross_section_O2, absorption_cross_section_O = np.loadtxt("Data\phot_abs.dat", skiprows=8, unpack=True)
wavelength_fism, irradiance = np.loadtxt("Data\Fism_daily_hr19990216.dat",  unpack = True, skiprows = 51, max_rows=949, delimiter=',', usecols = (1, 2))
electron_density, ion_temp, electron_temp = np.loadtxt("Data\IRI.dat", unpack = True, skiprows= 46, usecols=(1, 2, 3))


# finds T_r by taking the mean of ion_temp and nutral_temp
T_r = (ion_temp + nutral_temp)/2
print(T_r)

""" 
Task 0: integrate the continuity-equations for 3600 seconds with a constant ionization-rate of
        1*10^8 (/m^3/s) to obtain a stable background to use for the initial conditions of the different densities.
"""



""" 
Task 1: The response to an ionization pulse with 100 s duration should be modelled over a 600 s
        long time-period. The ion-electron production rate should be set to 1*10^10(m^3s^1) for 100 s
        followed by absolutely no further ionization. This modelling should be done at altitudes of 110,
        170 and 230 km.
"""