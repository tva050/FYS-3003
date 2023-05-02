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

# Altitude for the production rate


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
def production_rate(I, wavelength, threshold, sigma, n, alititude):
    h = 6.626e-34 # J s planck constant
    c = 3e8 # m/s speed of light
    eV = 6.242e18 # 1 eV = 6.242e18 J

    E = h*c * (1 / wavelength[wavelength <= threshold] - 1/threshold)*eV
    
    P = n[alititude] * sigma * I[alititude]
    P = P[wavelength <= threshold]*1e-6 # -> cm^-3 s^-1 
    return P, E

def phot_e_prod(altitude, all_I):
    P_N2, E_N2 = production_rate(all_I, wavelength_fism, N2_threshold, absorption_cross_section_N2, N2, altitude)
    P_O2, E_O2 = production_rate(all_I, wavelength_fism, O2_threshold, absorption_cross_section_O2, O2, altitude)
    P_O, E_O = production_rate(all_I, wavelength_fism, O_threshold, absorption_cross_section_O, O, altitude)
    
    energy = np.linspace(np.min(E_O), np.max(E_O2), 1000)

    P_N2_interpolated = sc.interpolate.interp1d(E_N2, P_N2, fill_value = "extrapolate", kind='linear')
    P_O2_interpolated = sc.interpolate.interp1d(E_O2, P_O2, fill_value = "extrapolate",kind='linear')
    P_O_interpolated  = sc.interpolate.interp1d(E_O, P_O, fill_value = "extrapolate",kind='linear')

    P_N2 = P_N2_interpolated(energy)
    P_O2 = P_O2_interpolated(energy)
    P_O  = P_O_interpolated(energy)

    p = P_N2 + P_O2 + P_O
    return p, energy 

""" energies = phot_e_prod(0, all_I[0])[1]
p = np.zeros((len(height), len(energies)))
for z in range(len(height)): 
    p[z] = phot_e_prod(z, all_I[0])[0]

plt.pcolormesh(energies, height*1e-3, p, norm = colors.LogNorm(vmin = 1e3, vmax = 1e-1))
plt.xlabel("Energy [eV]")
plt.ylabel("Altitude [km]")
plt.title("Production rate of photo-electrons SZA = 75$^\circ$")
plt.colorbar(label = "Production rate [cm$^{-3}$ s$^{-1}$]")
plt.show() """

p100_sza0 = phot_e_prod(150, all_I[0])[0]
energy100_sza0 = phot_e_prod(150, all_I[0])[1]
p100_sza20 = phot_e_prod(150, all_I[2])[0]
energy100_sza20 = phot_e_prod(150, all_I[2])[1]
p100_sza40 = phot_e_prod(150, all_I[4])[0]
energy100_sza40 = phot_e_prod(150, all_I[4])[1]
p100_sza75 = phot_e_prod(150, all_I[6])[0]
energy100_sza75 = phot_e_prod(150, all_I[6])[1]

p200_sza0 = phot_e_prod(200, all_I[0])[0]
energy200_sza0 = phot_e_prod(200, all_I[0])[1]
p200_sza20 = phot_e_prod(200, all_I[2])[0]
energy200_sza20 = phot_e_prod(200, all_I[2])[1]
p200_sza40 = phot_e_prod(200, all_I[4])[0]
energy200_sza40 = phot_e_prod(200, all_I[4])[1]
p200_sza75 = phot_e_prod(200, all_I[6])[0]
energy200_sza75 = phot_e_prod(200, all_I[6])[1]

p350_sza0 = phot_e_prod(350, all_I[0])[0]
energy350_sza0 = phot_e_prod(350, all_I[0])[1]
p350_sza20 = phot_e_prod(350, all_I[2])[0]
energy350_sza20 = phot_e_prod(350, all_I[2])[1]
p350_sza40 = phot_e_prod(350, all_I[4])[0]
energy350_sza40 = phot_e_prod(350, all_I[4])[1]
p350_sza75 = phot_e_prod(350, all_I[6])[0]
energy350_sza75 = phot_e_prod(350, all_I[6])[1]

plt.plot(energy100_sza0, p100_sza0, label = "SZA = 0$^\circ$")
plt.plot(energy100_sza20, p100_sza20, label = "SZA = 20$^\circ$", color = "gray")
plt.plot(energy100_sza40, p100_sza40, label = "SZA = 40$^\circ$", color = "orange")
plt.plot(energy100_sza75, p100_sza75, label = "SZA = 75$^\circ$", color = "teal")
plt.xlabel("Energy [eV]")
plt.ylabel("Production rate [cm$^{-3}$ s$^{-1}$]")
plt.xscale("log")
plt.yscale("log")
plt.title("Production rate of photo-electrons at 150 km altitude")
plt.legend()
plt.show()

plt.plot(energy200_sza0, p200_sza0, label = "SZA = 0$^\circ$")
plt.plot(energy200_sza20, p200_sza20, label = "SZA = 20$^\circ$", color = "gray")
plt.plot(energy200_sza40, p200_sza40, label = "SZA = 40$^\circ$", color = "orange")
plt.plot(energy200_sza75, p200_sza75, label = "SZA = 75$^\circ$", color = "teal")
plt.xlabel("Energy [eV]")
plt.ylabel("Production rate [cm$^{-3}$ s$^{-1}$]")
plt.xscale("log")
plt.yscale("log")
plt.title("Production rate of photo-electrons at 200 km altitude")
plt.legend()
plt.show()

plt.plot(energy350_sza0, p350_sza0, label = "SZA = 0$^\circ$")
plt.plot(energy350_sza20, p350_sza20, label = "SZA = 20$^\circ$", color = "gray")
plt.plot(energy350_sza40, p350_sza40, label = "SZA = 40$^\circ$", color = "orange")
plt.plot(energy350_sza75, p350_sza75, label = "SZA = 75$^\circ$", color = "teal")
plt.xlabel("Energy [eV]")
plt.ylabel("Production rate [cm$^{-3}$ s$^{-1}$]")
plt.xscale("log")
plt.yscale("log")
plt.title("Production rate of photo-electrons at 350 km altitude")
plt.legend()
plt.show()


""" altitudes = np.array([125, 140, 204, 392]) # Same altitudes as in figure 2.5.1 in M H Rees
energye = np.linspace(np.min(energies), np.max(energies), 1000)

plt.plot(energies, p[altitudes[0]], label = "Altitude = 125 km")
plt.plot(energies, p[altitudes[1]], label = "Altitude = 140 km", color = "gray")
plt.plot(energies, p[altitudes[2]], label = "Altitude = 204 km", color = "orange")
plt.plot(energies, p[altitudes[3]], label = "Altitude = 392 km", color = "teal")
plt.xlabel("Energy [eV]")
plt.ylabel("Production rate [cm$^{-3}$ s$^{-1}$]")
plt.title("Production rate of photo-electrons as a function of energy")
plt.xscale("log")
plt.yscale("log")
plt.legend() """
#plt.show()



# b) Make functions that calculate the photo-ionization profiles as a function of altitude (Eq. 2.3.4 M H Rees)
def ionization_profile(I, wavelength, threshold, sigma, n):
    q = np.zeros(height.shape[0]) 
    for i in range(height.shape[0]):
        wavelength_ = wavelength <= threshold
        integrand = I[i, wavelength_] * sigma[wavelength_]
        integral  = sc.integrate.trapz(integrand, wavelength[wavelength_])
        q[i] = integral * n[i]
    return q*1e9

q_N2 = ionization_profile(all_I[6], wavelength_fism, N2_threshold, absorption_cross_section_N2, N2)*1e-6
q_O2 = ionization_profile(all_I[6], wavelength_fism, O2_threshold, absorption_cross_section_O2, O2)*1e-6
q_O = ionization_profile(all_I[6], wavelength_fism, O_threshold, absorption_cross_section_O, O)*1e-6
total_q = q_N2 + q_O2 + q_O

plt.plot(q_O, height*1e-3, label="O")
plt.plot(q_O2, height*1e-3, label="O2", color = "orange")
plt.plot(q_N2, height*1e-3, label="N2", color = "gray")
plt.plot(total_q, height*1e-3, label="Total", color = "teal")
plt.xlabel("q [cm$^{-3}$ s$^{-1}$]")
plt.ylabel("Height [km]")
plt.title("photo-ionization rate SZA = 75$^\circ$")
plt.legend()
plt.show()


# Chapman 

def chapman(q, H, chi):
    q_m0 = np.max(q)
    z_m0 = height[np.argmax(q)]
    
    
