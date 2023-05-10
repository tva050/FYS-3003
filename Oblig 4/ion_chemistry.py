import numpy as np
import matplotlib.pyplot as plt

plt.style.use('ggplot')

height, O, N2, O2, mass_density, nutral_temp = np.loadtxt(r'Data\MSIS.dat', unpack = True ,skiprows=118, usecols= (0, 1, 2, 3, 4, 5))
electron_density, ion_temp, electron_temp, O_ions, H_ions, He_ions, O2_ions, No_ions, cluster_ions, N_ions  = np.loadtxt("Data\IRI.dat", unpack = True, skiprows= 46, usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))


# finds T_r by taking the mean of ion_temp and nutral_temp
T_r = (ion_temp + nutral_temp)/2
T_e = electron_temp

t = 3600 # seconds of integration


""" 
Task 0: integrate the continuity-equations for 3600 seconds with a constant ionization-rate of
        1*10^8 (/m^3/s) to obtain a stable background to use for the initial conditions of the different densities.
"""
n_O_ions = O_ions /       # calculate the number density of 0 which is given with percent

def myODE(x, t, altitudes): 
        
        # constants ish
        alpha1 = 2.1e-13 * (T_e/300)**-0.85
        alpha2 = 1.9e-13 * (T_e/300)**-0.5
        alpha3 = 1.8e-13 * (T_e/300)**-0.39
        
        alpha_r = 3.7e-18 * (250/T_e)**-0.7
        
        k1 = 2e-18
        k2 = 2e-17 * (T_r/300)**-0.4
        k3 = 4.4e-16
        k4 = 5e-22
        k5 = 1.4e-16 * (T_r/300)**-0.44
        k6 = 5e-17 * (T_r/300)**-0.8
        
        
        #def solver(x0, t):
                



""" 
Task 1: The response to an ionization pulse with 100 s duration should be modelled over a 600 s
        long time-period. The ion-electron production rate should be set to 1*10^10(m^3s^1) for 100 s
        followed by absolutely no further ionization. This modelling should be done at altitudes of 110,
        170 and 230 km.
"""


""" 
Task 2: Same as above but increase the electron and ion temperatures by 1000 K below 150 km,
        and by 2000 K at the altitudes above 150 km. Compare the ion-composition for the two cases.
"""