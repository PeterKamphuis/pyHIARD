#-*- coding: future_fstrings -*-


global H_0
H_0 = 69.7  # km/s/Mpc
global c
c = 299792458  # light speed in m/s
global c_kms
c_kms = 299792.458  # light speed in km/s
global pc
pc = 3.086e+18  # parsec in cm
global solar_mass
solar_mass = 1.98855e30  # Solar mass in kg
global solar_luminosity
solar_luminosity = 3.828e26    # Bolometric Solar Luminosity in W
global HI_mass
HI_mass = 1.6737236e-27  # Mass of hydrogen in kg
global G
G = 6.67430e-11  # m^3/kg*s^2
global Gsol
Gsol = G/(1000.**3)*solar_mass  # km^3/M_sol*s^2
# transform (pc x km^2)/(s^2 x solarmass)
global G_agc
G_agc = Gsol/(pc/(1000.*100.))
global HI_rest_freq
HI_rest_freq=1.4204057517667e+9  # Hz

# G=6.674 x 10^neg20 #km^3⋅kg^neg1⋅s^neg2 pc=3.086e+13 km solarmass=1.98855e30 kg
# transform (pc x km^2)/(s^2 x solarmass)
#G= 4.30058*10**-3
