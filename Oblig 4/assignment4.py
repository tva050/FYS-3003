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

""" 
Task 0: integrate the continuity-equations for 3600 seconds with a constant ionization-rate of
        1*10^8 (/m^3/s) to obtain a stable background to use for the initial conditions of the different densities.
"""

def task_0(altitude):
    # finds T_r by taking the mean of ion_temp and nutral_temp
    T_r = (ion_temp + nutral_temp)/2
    T_e = electron_temp
    
    def ODEs(ni, t, altitude):
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
    
        # Ionization rates for the different species
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
        d_nNO = -k3*O2p*No + k4*O2p*nN2
        d_ne = q_e - alpha1*NOp*ne - alpha2*O2p*ne - alpha3*N2p*ne - alpha_r*Op*ne
    
        if not np.allclose(q_e, q_e): # if q_e is not constant raise an error
            raise ValueError("q_e is not constant")
    
        return [d_nN2p[altitude], d_nO2p[altitude], d_nOp[altitude], d_nNOp[altitude], d_nNO[altitude], d_ne[altitude]]

    # Initial conditions
        # nN2p, nO2p,         nOp,         nNOp,      nNO, ne
    ni0 = [0 , n_O2_ions[0], n_O_ions[0], n_NO_ions[0], 0, n_e[0]]
    t = np.linspace(0, 3600, 3600) # time vector 
    sol1 = odeint(ODEs, ni0, t, args= (altitude,)) # solving the ODEs for the wanted altitude

    # Extracting the different densities from the solution
    N2p = sol1[:,0] 
    O2p = sol1[:,1]
    Op  = sol1[:,2]
    NOp = sol1[:,3]
    NO  = sol1[:,4]
    ne  = sol1[:,5]

    plt.plot(t, N2p, label = "$N_2^+$")
    plt.plot(t, O2p, label = "$O_2^+$")
    plt.plot(t, Op, label = "$O^+$")
    plt.plot(t, NOp, label = "$NO^+$")
    plt.plot(t, NO, label = "$NO$")
    plt.plot(t, ne, label = "$e^- $")
    plt.xlabel("Time [s]")
    plt.ylabel("Density [m^-3]")
    plt.yscale("log")
    plt.title("Ion densities at " + str(int(altitude)) + "km, $q_e = 1\cdot 10^8$")
    plt.legend()
    plt.show()

        
""" 
Task 1: The response to an ionization pulse with 100 s duration should be modelled over a 600 s
        long time-period. The ion-electron production rate should be set to 1*10^10(m^3s^1) for 100 s
        followed by absolutely no further ionization. This modelling should be done at altitudes of 110,
        170 and 230 km.
"""

def task_1(altitude):
    # finds T_r by taking the mean of ion_temp and nutral_temp
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
    time = np.linspace(0, 600, 600) # time vector
    solve = odeint(ODEs, ni0, time, args= (altitude,)) # solving the ODEs for the wanted altitude

    N2p = solve[:,0] 
    O2p = solve[:,1]
    Op  = solve[:,2]
    NOp = solve[:,3]
    NO  = solve[:,4]
    ne  = solve[:,5]

    plt.plot(time, N2p, label = "$N_2^+$")
    plt.plot(time, O2p, label = "$O_2^+$")
    plt.plot(time, Op, label = "$O^+$")
    plt.plot(time, NOp, label = "$NO^+$")
    plt.plot(time, NO, label = "$NO$")
    plt.plot(time, ne, label = "$e^- $")
    plt.xlabel("Time [s]")
    plt.ylabel("Density [m^-3]")
    plt.yscale("log")
    plt.ylim(1e7, 1e12)
    plt.title("Ion densities at " + str(int(altitude)) + "km, $q_e = 1\cdot 10^{10}$")
    plt.legend()
    plt.show()

""" 
Task 2: Same as above but increase the electron and ion temperatures by 1000 K below 150 km,
        and by 2000 K at the altitudes above 150 km. Compare the ion-composition for the two cases.
"""
def task_2(altitude):
    # finds T_r by taking the mean of ion_temp and nutral_temp
    T_e = electron_temp
    
    # function to increase the temperature
    def temperature_increase(altitude):
        # 
        if altitude < 150:
            return T_e + 1000, ion_temp + 1000
        else: 
            return T_e + 2000, ion_temp + 2000
    
    T_e = temperature_increase(altitude)[0]
    T_r = (temperature_increase(altitude)[1] + nutral_temp)/2
    
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
    time = np.linspace(0, 600, 600) # time vector
    solve = odeint(ODEs, ni0, time, args= (altitude,)) # solving the ODEs for the wanted altitude

    N2p = solve[:,0] 
    O2p = solve[:,1]
    Op  = solve[:,2]
    NOp = solve[:,3]
    NO  = solve[:,4]
    ne  = solve[:,5]

    plt.plot(time, N2p, label = "$N_2^+$")
    plt.plot(time, O2p, label = "$O_2^+$")
    plt.plot(time, Op, label = "$O^+$")
    plt.plot(time, NOp, label = "$NO^+$")
    plt.plot(time, NO, label = "$NO$")
    plt.plot(time, ne, label = "$e^- $")
    plt.xlabel("Time [s]")
    plt.ylabel("Density [m^-3]")
    plt.yscale("log")
    plt.ylim(1e7, 1e12)
    plt.title("Ion densities at " + str(int(altitude)) + "km, $q_e = 1\cdot 10^{10}$")
    plt.legend()
    plt.show()
    
""" 
Task 3: Compare the electron-density decay at 110 and 230 km with the expected decrease-characteristics:
        "alpha" ((ne(t_off )/(1 + alpha_e*ne(t_off ))) and "beta" (exponential decay). Remember to
        include both the dissociative recombination of O2p and NOp when calculating the beta-decay.
"""
def task_3(altitude): 
    T_r = (ion_temp + nutral_temp)/2
    T_e = electron_temp
    
    def ODEs(ni, t, altitude):
        global alpha1, alpha2, alpha3, k1, k2
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
    
        # Ionization rates for the different species
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
        d_nNO = -k3*O2p*No + k4*O2p*nN2
        d_ne = q_e - alpha1*NOp*ne - alpha2*O2p*ne - alpha3*N2p*ne - alpha_r*Op*ne
    
        if not np.allclose(q_e, q_e): # if q_e is not constant raise an error
            raise ValueError("q_e is not constant")
    
        return [d_nN2p[altitude], d_nO2p[altitude], d_nOp[altitude], d_nNOp[altitude], d_nNO[altitude], d_ne[altitude]]

    # Initial conditions
        # nN2p, nO2p,         nOp,         nNOp,      nNO, ne
    ni0 = [0 , n_O2_ions[0], n_O_ions[0], n_NO_ions[0], 0, n_e[0]]
    t = np.linspace(0, 3600, 3600) # time vector 
    sol1 = odeint(ODEs, ni0, t, args= (altitude,)) # solving the ODEs for the wanted altitude

    # Extracting the different densities from the solution
    N2p = sol1[:,0] 
    O2p = sol1[:,1]
    Op  = sol1[:,2]
    NOp = sol1[:,3]
    NO  = sol1[:,4]
    ne  = sol1[:,5]

    t_off = 0 
    
    alpha_e = alpha1*(NOp/ne) + alpha2*(O2p/ne) + alpha3*(N2p/ne) # refered to as alpha overline 
    def alpha(t_off):
        alpha = (ne[t_off])/(1 + alpha_e*ne[t_off])
        return alpha
    
    beta_marked = (k1*nN2 + k2*nO2) / (1 + (k1/alpha1)*(nN2/ne) + (k2/alpha2)*(nO2/ne))
    def beta_decay(t_off):
        beta_decay = ne[t_off]*np.exp(-beta_marked*(t-t_off))
        return beta_decay
    
    plt.plot(t, ne, label = "$n_e$")
    plt.plot(t, alpha(t_off), label = "$\\alpha$")
    plt.plot(t, beta_decay(t_off), label = "$\\beta$")
    plt.xlabel("Time [s]")
    plt.ylabel("Density [m^-3]")
    plt.yscale("log")
    plt.title("Electron density at " + str(int(altitude)) + "km")
    plt.legend()
    plt.show()
    
    
    
if __name__ == "__main__":
    # !!! Remember to run the functions, the altitudes need to be indexed !!!
    # The altitudes to choose from are 110 ([0]), 170 ([1]) and 230([2]) km
    
    #task_0(altitudes[0])
    #task_1(altitudes[1])
    #task_2(altitudes[1])
    task_3(altitudes[0])
