import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

plt.style.use('ggplot')

""" Constants """
k_B = 1.38064852e-23 # Boltzmann constant (J/K)
R_E = 6_371e6 # Radius of Earth (m)
g_0 = 9.80665 # Gravitational acceleration at sea level (m/s2)
P_0 = 101_325 # Pressure at sea level (Pa)
# Atmoic mass of O, N2, O2 (kg)
m_O = 15.9994e-3
m_N2 = 28.0134e-3
m_O2 = 31.9988e-3

 
""" Task 1 """

# Read data from file MSIS.dat in a reasonable format
# height (km), O (cm-3), N2 (cm-3), O2 (cm-3), mass_density (g/cm3), temperature_neutral (K)
height, O, N2, O2, mass_density, temperature_neutral = np.loadtxt(r'Data\MSIS.dat', skiprows=18, unpack=True)

# total number density (cm-3) convert to m-3
number_density = (O + N2 + O2)*1e6
# mass density (kg/m3)
mass_density = mass_density * 1e3 

height = height * 1e3 # convert to m

#ave_mass = (O * m_O + N2 * m_N2 + O2 * m_O2) / (O + N2 + O2) 
ave_mass = mass_density / number_density
# Gravitational accereleration which varies with height.
def gravitational_acceleration(height):
    g = g_0 *(R_E**2 / (R_E + height)**2)
    return g

# Atmosphere-alt-vars.pdf page 5, for the scale height (H) equation.
def scale_height(temperature, ave_mass):
    return k_B * temperature / (ave_mass * gravitational_acceleration(height))


plt.plot(scale_height(temperature_neutral, ave_mass)*1e-3, height*1e-3 ) # convert to km
plt.xlabel('Scale height (km)')
plt.ylabel('Height (km)')
plt.title('Scale height variation with height')
plt.show()


""" Task 2 """
# Using the barometric equation, to calculate the pressure variation with height.
def pressure_variation(height, pressure_0):
    pressure = pressure_0 * np.exp(-height / scale_height(temperature_neutral, ave_mass)*1e-3) # convert to km
    return pressure
plt.plot(pressure_variation(height, P_0)*1e-3, height*1e-3) 
plt.xlabel('Pressure (kPa)')
plt.ylabel('Height (km)')
plt.title('Pressure variation with height')
plt.show()




""" Task 3 """

# Ideal gas law, p(z) = n(z)k_BT(z)
# n(z) = n(z_0)*(T(z_0)/T(z))*exp(-z/H(z))

# Density variation with height
def density(height):
    density = number_density[0]* (temperature_neutral[0]/temperature_neutral) * np.exp(-height/scale_height(temperature_neutral, ave_mass)*1e-3)
    return density
    
def ideal_gaslaw(height):
    ideal_gas_law = density(height) * k_B * temperature_neutral
    return ideal_gas_law

plt.subplot(1,2,1)
plt.plot(pressure_variation(height, P_0)*1e-3, height*1e-3)
plt.xlabel('Pressure (kPa)')
plt.ylabel('Height (km)')
plt.title('Pressure variation')

plt.subplot(1,2,2)
plt.plot(ideal_gaslaw(height)*1e-3, height*1e-3)
plt.xlabel('Pressure (kPa)')
plt.title('Ideal gas law')

plt.show()


def temperature_vs_adiabatic():
    dT = np.diff(temperature_neutral) # Temperature difference
    dz = np.diff(height) # Height difference
    
    dTdz = dT/dz # Temperature gradient
    
    abc_lapsrate = -(((1.4-1)/(1.4)) * ((ave_mass *g_0)/k_B))*1e3 # Adiabatic lapse rate, with coneversion to K/km
    
    plt.plot(height[:-1], dTdz)
    plt.plot(height, abc_lapsrate)
    plt.xlabel('Height (km)')
    plt.ylabel('Temperature gradient (K/km)')
    plt.show()
    
temperature_vs_adiabatic()
    
