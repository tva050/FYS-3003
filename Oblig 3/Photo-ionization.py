import matplotlib.pyplot as plt
from matplotlib import colors
import os
import scipy as sc
import numpy as np

""" ---------------------------- Import data ---------------------------- """
I_SZA0  = np.loadtxt(os.path.join("Data", "I_SZA0.csv"), unpack = True, delimiter=",")
I_SZA10 = np.loadtxt(os.path.join("Data", "I_SZA10.csv"), unpack = True, delimiter=",")
I_SZA20 = np.loadtxt(os.path.join("Data", "I_SZA20.csv"), unpack = True, delimiter=",")
I_SZA30 = np.loadtxt(os.path.join("Data", "I_SZA30.csv"), unpack = True, delimiter=",")
I_SZA40 = np.loadtxt(os.path.join("Data", "I_SZA40.csv"), unpack = True, delimiter=",")
I_SZA50 = np.loadtxt(os.path.join("Data", "I_SZA50.csv"), unpack = True, delimiter=",")
I_SZA75 = np.loadtxt(os.path.join("Data", "I_SZA75.csv"), unpack = True, delimiter=",")

all_I = np.array([I_SZA0, I_SZA10, I_SZA20, I_SZA30, I_SZA40, I_SZA50, I_SZA75])
szas = np.array([0, 10, 20, 30, 40, 50, 75])

height, O, N2, O2 = np.loadtxt(r'Data\MSIS.dat', unpack = True ,skiprows=18, usecols= (0, 1, 2, 3))
wavelength, absorption_cross_section_N2, absorption_cross_section_O2, absorption_cross_section_O = np.loadtxt("Data\phot_abs.dat", skiprows=8, unpack=True)
wavelength_fism, irradiance = np.loadtxt("Data\Fism_daily_hr19990216.dat",  unpack = True, skiprows = 51, max_rows=949, delimiter=',', usecols = (1, 2))

# Converting to SI units
wavelength_fism *= 1e-9 # m
height *= 1e3

O, N2, O2 = O/1e6, N2/1e6, O2/1e6

# threshold wavelengths
wavelength_treshold = 0.1e-6 # m

# ionization threshold
N2_treshold = 7.96e-8 # Ångström -> m
O2_treshold = 1.026e-7 # Ångström -> m
O_treshold = 9.11e-8 # Ångström -> m



# a) Make functions that calculate the production rate of photo-electrons as a function of altitude and energy.


# b) Make functions that calculate the photo-ionization profiles as a function of altitude

def ionisation_profile(I, wavelength, treshold, sigma, N):
    q = np.zeros(height.shape[0])
    for i in range(height.shape[0]):
        integrand = I[i, wavelength <= treshold] * sigma[wavelength <= treshold]
        integral  = np.trapz(integrand, wavelength[wavelength <= treshold])
        q[i] = integral * N[i]
    return q*1e9

q_O = ionisation_profile(all_I[0], wavelength_fism, O_treshold, absorption_cross_section_O, O)*1e-6

plt.plot(q_O, height*1e-3)
plt.xlabel("q [cm$^{-3}$ s$^{-1}$]")
plt.ylabel("Height [km]")
plt.title("O ionisation profile")
plt.show()



