import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

plt.style.use('ggplot')

height, nO, nN2, nO2, mass_density, nutral_temp = np.loadtxt(r'Data\MSIS.dat', unpack = True ,skiprows=118, usecols= (0, 1, 2, 3, 4, 5))
n_e, ion_temp, electron_temp, O_ions, H_ions, He_ions, O2_ions, No_ions, cluster_ions, N_ions  = np.loadtxt("Data\IRI.dat", unpack = True, skiprows= 46, usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

nO, nN2, nO2 = nO*1e6, nN2*1e6, nO2*1e6 # converting to /m^3

# finds T_r by taking the mean of ion_temp and nutral_temp
T_r = (ion_temp + nutral_temp)/2
T_e = electron_temp

""" 
Task 0: integrate the continuity-equations for 3600 seconds with a constant ionization-rate of
        1*10^8 (/m^3/s) to obtain a stable background to use for the initial conditions of the different densities.
"""

# constants ish
sum_ions = O_ions + O2_ions + No_ions + N_ions
n_O_ions = (O_ions/sum_ions)*n_e
n_O2_ions = (O2_ions/sum_ions)*n_e
n_NO_ions = (No_ions/sum_ions)*n_e
n_N_ions = (N_ions/sum_ions)*n_e
other_ions = H_ions + He_ions + cluster_ions

def ODEs(ni, t, altitude):  
    
    
    alpha1 = 2.1e-13 * (T_e/300)**-0.85
    alpha2 = 1.9e-13 * (T_e/300)**-0.5
    alpha3 = 1.8e-13 * (T_e/300)**-0.39

    alpha_r = 3.7e-18 * (250/T_e)**0.7

    k1 = 2e-18
    k2 = 2e-17 * (T_r/300)**-0.4
    k3 = 4.4e-16
    k4 = 5e-22
    k5 = 1.4e-16 * (T_r/300)**-0.44
    k6 = 5e-17 * (T_r/300)**-0.8

    q_e =  1e8
    q_N2 = q_e * (0.92 * nN2) / (0.92*nN2 + nO2+ 0.56 * nO)
    q_O2 = q_e * (nO2) / (0.92*nN2 + nO2+ 0.56 * nO)
    q_O =  q_e * (0.56 * nO) / (0.92*nN2 + nO2+ 0.56 * nO)
    
    #------------------------------------------- ODES -------------------------------------------#
    
    # Assign each ODE to a vector element 
    N2p = ni[0]
    O2p = ni[1]
    Op = ni[2]
    NOp = ni[3]
    No = ni[4]
    ne = ni[5]
    
    d_nN2p = q_N2 - alpha3*N2p*ne - k5*N2p*nO - k6*N2p*nO2
    d_nO2p = q_O2 - alpha2*O2p*ne + k2*nO2*Op - k3*O2p*No - k4*O2p*nN2 + k6*nO2*N2p
    d_nOp = q_O - k1*Op*nN2 - k2*Op*nO2 - alpha_r*Op*ne
    d_nNOp = -alpha1*NOp*ne + k1*Op*nN2 + k3*O2p*No + k4*O2p*nN2 + k5*N2p*nO
    d_nNO = -k3*O2p*No + k4*O2p*nN2 + k5*N2p*nO
    d_ne = q_e - alpha1*NOp*ne - alpha2*O2p*ne - alpha3*N2p*ne - alpha_r*Op*ne
    
    return [d_nN2p, d_nO2p, d_nOp, d_nNOp, d_nNO, d_ne]




# Initial conditions
    # nN2p, nO2p,         nOp,     nNOp,  nNO, ne
ni0 = [0 , n_O2_ions[0], n_O_ions[0], n_NO_ions[0], 0, n_e[0]]
altitudes = [110, 170, 230] # wanted altitudes in km
t = np.linspace(0, 3600, 3600) # time vector
print(t) 
sol = odeint(ODEs, ni0, t, args= (110,)) # solving the ODEs  


A = sol[:,0]
B = sol[:,1]
C = sol[:,2]
D = sol[:,3]
E = sol[:,4]

plt.plot(t, A, label = "N2+")
plt.plot(t, B, label = "O2+")
plt.plot(t, C, label = "O+")
plt.plot(t, D, label = "NO+")
plt.plot(t, E, label = "NO")
plt.xlabel("Time [s]")
plt.ylabel("Density [m^-3]")
plt.title("Ion densities at 110 km")
plt.legend()
plt.show()
