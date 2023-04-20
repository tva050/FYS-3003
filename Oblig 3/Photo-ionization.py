import matplotlib.pyplot as plt
from matplotlib import colors
import os
import scipy as sc
import numpy as np

plt.style.use('ggplot')

""" ---------------------------- Import data ---------------------------- """
I_SZA0  = np.loadtxt(os.path.join("Data", "I_SZA0.csv"), delimiter=",")
I_SZA10 = np.loadtxt(os.path.join("Data", "I_SZA10.csv"), delimiter=",")
I_SZA20 = np.loadtxt(os.path.join("Data", "I_SZA20.csv"), delimiter=",")
I_SZA30 = np.loadtxt(os.path.join("Data", "I_SZA30.csv"), delimiter=",")
I_SZA40 = np.loadtxt(os.path.join("Data", "I_SZA40.csv"), delimiter=",")
I_SZA50 = np.loadtxt(os.path.join("Data", "I_SZA50.csv"), delimiter=",")
I_SZA75 = np.loadtxt(os.path.join("Data", "I_SZA75.csv"), delimiter=",")

all_I = np.array([I_SZA0, I_SZA10, I_SZA20, I_SZA30, I_SZA40, I_SZA50, I_SZA75])*1e4
szas = np.array([0, 10, 20, 30, 40, 50, 75])


height, O, N2, O2 = np.loadtxt(r'Data\MSIS.dat', unpack = True ,skiprows=18, usecols= (0, 1, 2, 3))
wavelength, absorption_cross_section_N2, absorption_cross_section_O2, absorption_cross_section_O = np.loadtxt("Data\phot_abs.dat", skiprows=8, unpack=True)
wavelength_fism, irradiance = np.loadtxt("Data\Fism_daily_hr19990216.dat",  unpack = True, skiprows = 51, max_rows=949, delimiter=',', usecols = (1, 2))

# Converting to SI units
wavelength_fism *= 1e-9 # m
height *= 1e3

O, N2, O2 = O*1e6, N2*1e6, O2*1e6

# threshold wavelengths
wavelength_treshold = 0.1e-6 # m

# ionization threshold
N2_threshold = 7.96e-8 # Ångström -> m
O2_threshold = 1.026e-7 # Ångström -> m
O_threshold = 9.11e-8 # Ångström -> m


#___________________________________________________________
#Since data have different sizes, we need to interpolate the data to the same size
interpolate_N2 = sc.interpolate.interp1d(wavelength, absorption_cross_section_N2, kind='linear')
interpolate_O2 = sc.interpolate.interp1d(wavelength, absorption_cross_section_O2, kind='linear')
interpolate_O  = sc.interpolate.interp1d(wavelength, absorption_cross_section_O,  kind='linear')

absorption_cross_section_N2 = interpolate_N2(wavelength_fism)
absorption_cross_section_O2 = interpolate_O2(wavelength_fism)
absorption_cross_section_O  = interpolate_O(wavelength_fism)
#___________________________________________________________


# a) Make functions that calculate the production rate of photo-electrons as a function of altitude and energy.
def production_rate(I, wavelength, threshold, sigma, n):
    h = 6.626e-34 # J s planck constant
    c = 3e8 # m/s speed of light
    eV = 6.242e18 # 1 eV = 6.242e18 J

    E = h*c * (1 / wavelength[wavelength <= threshold] - 1/threshold)*eV
    
    P = n[200] * sigma * I[200]
    P = P[wavelength <= threshold]*1e-6
    return P, E

P_N2, E_N2 = production_rate(all_I[0], wavelength_fism, N2_threshold, absorption_cross_section_N2, N2)
P_O2, E_O2 = production_rate(all_I[0], wavelength_fism, O2_threshold, absorption_cross_section_O2, O2)
P_O, E_O = production_rate(all_I[0], wavelength_fism, O_threshold, absorption_cross_section_O, O)

energy = np.linspace(np.min(E_O), np.max(E_O2), 1000)

P_N2_interpolated = sc.interpolate.interp1d(E_N2, P_N2, fill_value =  ,kind='linear')
P_O2_interpolated = sc.interpolate.interp1d(E_O2, P_O2, kind='linear')
P_O_interpolated  = sc.interpolate.interp1d(E_O, P_O, kind='linear')

P_N2 = P_N2_interpolated(energy)
P_O2 = P_O2_interpolated(energy)
P_O  = P_O_interpolated(energy)

plt.plot(energy, P_N2 + P_O2 + P_O)
plt.show()

# b) Make functions that calculate the photo-ionization profiles as a function of altitude (Eq. 2.3.4 M H Rees)
def ionization_profile(I, wavelength, threshold, sigma, n):
    q = np.zeros(height.shape[0])
    for i in range(height.shape[0]):
        wavelength_ = wavelength <= threshold
        integrand = I[i, wavelength_] * sigma[wavelength_]
        integral  = sc.integrate.trapz(integrand, wavelength[wavelength_])
        q[i] = integral * n[i]
    return q*1e9

q_N2 = ionization_profile(all_I[0], wavelength_fism, N2_threshold, absorption_cross_section_N2, N2)*1e-6
q_O2 = ionization_profile(all_I[0], wavelength_fism, O2_threshold, absorption_cross_section_O2, O2)*1e-6
q_O = ionization_profile(all_I[0], wavelength_fism, O_threshold, absorption_cross_section_O, O)*1e-6
total_q = q_N2 + q_O2 + q_O

plt.plot(q_O, height*1e-3, label="O")
plt.plot(q_O2, height*1e-3, label="O2")
plt.plot(q_N2, height*1e-3, label="N2")
plt.plot(total_q, height*1e-3, label="Total", color = "teal")
plt.xlabel("q [cm$^{-3}$ s$^{-1}$]")
plt.ylabel("Height [km]")
plt.title("photo-ionization rate")
plt.legend()
plt.show()



