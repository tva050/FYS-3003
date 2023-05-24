import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

plt.style.use('ggplot')

height, nO, nN2, nO2, mass_density, nutral_temp = np.loadtxt(r'Data\MSIS.dat', unpack = True ,skiprows=118, usecols= (0, 1, 2, 3, 4, 5))
n_e, ion_temp, electron_temp, O_ions, H_ions, He_ions, O2_ions, No_ions, cluster_ions, N_ions  = np.loadtxt("Data\IRI.dat", unpack = True, skiprows= 46, usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

nO, nN2, nO2 = nO*1e6, nN2*1e6, nO2*1e6

altitudes = [110, 170, 230]

# constants ish
sum_ions = O_ions + O2_ions + No_ions + N_ions
n_O_ions = (O_ions/sum_ions)*n_e
n_O2_ions = (O2_ions/sum_ions)*n_e
n_NO_ions = (No_ions/sum_ions)*n_e
n_N_ions = (N_ions/sum_ions)*n_e
other_ions = H_ions + He_ions + cluster_ions

def task_3(altitude): 
    T_r = (ion_temp + nutral_temp)/2
    T_e = electron_temp
    
    # constants ish
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
        
    def ionization_rate(t):
        q_e = 1e10 # ionization rate
        if t < 100:
            return q_e
        else:
            return 0

    
    def ODEs(ni, t, altitude):
        #------------------------------------------- ODES -------------------------------------------#
    
        q_e = ionization_rate(t)
        q_N2 = q_e * (0.92 * nN2) / (0.92*nN2 + nO2+ 0.56 * nO)
        q_O2 = q_e * (nO2) / (0.92*nN2 + nO2+ 0.56 * nO)
        q_O =  q_e * (0.56 * nO) / (0.92*nN2 + nO2+ 0.56 * nO)
    
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
        d_nNO = -k3*O2p*No + k4*O2p*nN2 
        d_ne = q_e - alpha1*NOp*ne - alpha2*O2p*ne - alpha3*N2p*ne - alpha_r*Op*ne
    
        return [d_nN2p[altitude], d_nO2p[altitude], d_nOp[altitude], d_nNOp[altitude], d_nNO[altitude], d_ne[altitude]]

    # Initial conditions
        # nN2p, nO2p,         nOp,         nNOp,      nNO, ne
    ni0 = [0 , n_O2_ions[0], n_O_ions[0], n_NO_ions[0], 0, n_e[0]]
    t = np.linspace(0, 600, 601) # time vector
    solve = odeint(ODEs, ni0, t, args= (altitude,)) # solving the ODEs for the wanted altitude

    N2p = solve[:,0] 
    O2p = solve[:,1]
    Op  = solve[:,2]
    NOp = solve[:,3]
    NO  = solve[:,4]
    ne  = solve[:,5]
    
    
    beta = (k1*nN2[altitude] + k2[altitude]*nO2[altitude]) / (1 + (k1/alpha1[altitude])*(nN2[altitude]/ ne[100]) + (k2[altitude]/alpha2[altitude])*(nO2[altitude]/ ne[100]))
    alpha_eff = ((alpha1[altitude]/(k1*nN2[altitude])) + alpha2[altitude]/(k2[altitude]*nO2[altitude])) * beta
    
    #alpha_eff = alpha1[altitude]*(NOp[100]/ne[100]) + alpha2[altitude]*(O2p[100]/ne[100]) + alpha3[altitude]*(N2p[100]/ne[100])
    
    def beta_decay(t):
        return ne[100]*np.exp(-beta*t[100:])
        
    def alpha_decay(t):
        return ne[100]/(1 + alpha_eff*ne[100]*t[100:])



    plt.plot(t, ne, label = "$n_{e^-}$")
    plt.plot(t[100:], beta_decay(t), label = "$\\beta$ decay")
    plt.plot(t[100:], alpha_decay(t), label = "$\\alpha$ decay")
    plt.xlabel("Time [s]")
    plt.ylabel("Density [m^-3]")
    plt.yscale("log")
    #plt.ylim(2.9e10, 1e12)
    plt.title("Electron density at " + str(int(altitude)) + "km, $q_e = 1\cdot 10^{10}$")
    plt.legend()
    plt.show()
    
task_3(110)

def noe(var1, var2=100, var3=20):
    return print("var1  is", var1/2, "var2 is", var2, "and var3", var3)

