import numpy as np
from astropy.constants import M_sun, G, L_sun, M_earth
import astropy.units as u

if __name__ == '__main__':


    # mass of TRAPPIST-1 e
    M_e = 0.692 * M_earth

    # mass of TRAPPIST-1
    M_star = 0.089 * M_sun

        # Luminosity of TRAPPIST-1
    L = 0.00055 * L_sun


    # current TRAPPIST-1 e orbital parameters
    T_E = 6.101013 * u.day
    r_E = (G * (M_star + M_e) * T_E**2 / (4 * np.pi**2))**(1/3)
    
    print(f"TRAPPIST-1 e orbit in AU: {r_E.to(u.au):.5f}")



    # ORBITAL PERIOD and RADIUS CHANGE ----------------------------

    # The coherence time for TRAPPIST-1 e's orbital period should be 
    # equal to the travel time for the laser beam to be randomized
    # at the time of its arrival to maximize chances on avoiding a hit

    print("\n ------------------------------------------")
    print("Orbital period and radius change:")
    # coherence time = travel time
    tau = 42. * u.year

    # uncertainty we need to introduce in TRAPPIST-1 e's orbital period
    # to avoid a hit
    # delta_T = (T_E + delta_T)**2 / tau

    # transform to a quadratic equation

    # delta_T**2 + 2 * T_E * delta_T + T_E**2 - tau * delta_T = 0

    # delta_T**2 + (2 * T_E - tau) * delta_T + T_E**2 = 0

    # solve quadratic equation
    delta_T = (-2 * T_E + tau - np.sqrt((2 * T_E - tau)**2 - 4 * T_E**2)) / 2

    print(f"Desired uncertainty in TRAPPIST-1 e's orbital radius: {delta_T.to(u.day):.5f}")

    # desired TRAPPIST-1 e orbital parameters
    T_Enew = T_E + delta_T

    # apply Kepler's third law
    r_Enew = (G * (M_star + M_e) * T_Enew**2 / (4 * np.pi**2))**(1/3)

    # potential energy of a circular orbit
    U = - G * M_star * M_e / r_E

    print(f"Potential energy of TRAPPIST-1 e orbit: {U.to(u.Joule):.1e}")

    print(f"Desired TRAPPIST-1 e orbit in AU: {r_Enew.to(u.au):.5f}")

    # calculate the energy difference, i.e. the difference in potential energy
    E = G * M_star * M_e * ( 1/r_E - 1/r_Enew)

    print(f"Energy difference: {E.to(u.Joule):.1e}")

    # Caclute the power needed to provide this energy difference in 42 years
    P = (E / tau).to(u.Watt)

    print(f"Power needed: {P:.1e}")

    # Power in units of TRAPPIST-1 luminosity
    P_TRAPPIST = P / L

    print(f"Power needed in units of TRAPPIST-1 luminosity: {P_TRAPPIST*100:.2f} %")

    print("\n ------------------------------------------")
    print("Orbital inclination change:")


    # INCLINATION CHANGE ------------------------------------------

    # Now what power would be needed to change the inclination of our orbit by 1 deg
    delta_i = 1 * u.deg

    # TRAPPIST-1 e's current speed in its orbit at apohelion, where the inclination 
    # change would be most efficient
    V = (2 * np.pi * r_E / T_E).to(u.km/u.s)

    # delta V equation from Wikipedia, assuming eccentricity is zero
    delta_V = 2 * V * np.sin(delta_i / 2)

    print(f"Delta V needed: {delta_V:.1f}")

    # calculate the energy difference, i.e. the difference in kinetic energy
    E = 0.5 * M_e * delta_V**2

    print(f"Energy difference: {E.to(u.Joule):.1e}")

    # Caclute the power needed to provide this energy difference in 42 years, that is 84 thrusts of 1h each
    P = (E / (84 * u.h)).to(u.Watt)

    print(f"Power needed for a {delta_i} inclination change: {P:.1e}")


    # Power in units of TRAPPIST-1 luminosity
    P_TRAPPIST = P / L

    print(f"Power needed for a {delta_i} inclination change in units of TRAPPIST-1 luminosity: {P_TRAPPIST*100:.2f} %")


