import numpy as np
from astropy.constants import M_sun, M_earth, G
import astropy.units as u

if __name__ == '__main__':



    # ORBITAL PERIOD and RADIUS CHANGE ----------------------------

    # The coherence time for Earth's orbital period should be 
    # equal to the travel time for the laser beam to be randomized
    # at the time of its arrival to maximize chances on avoiding a hit


    # current Earth orbital parameters
    r_E = 1 * u.au
    T_E = 1 * u.year

    # coherence time = travel time
    tau = 42. * u.year

    # uncertainty we need to introduce in Earth's orbital period
    # to avoid a hit
    # delta_T = (T_E + delta_T)**2 / tau

    # transform to a quadratic equation

    # delta_T**2 + 2 * T_E * delta_T + T_E**2 - tau * delta_T = 0

    # delta_T**2 + (2 * T_E - tau) * delta_T + T_E**2 = 0

    # solve quadratic equation
    delta_T = (-2 * T_E + tau - np.sqrt((2 * T_E - tau)**2 - 4 * T_E**2)) / 2

    print(f"Desired uncertainty in Earth's orbital radius: {delta_T.to(u.day):.1f}")

    # desired Earth orbital parameters
    T_Enew = 1 * u.year + delta_T

    # apply Kepler's third law
    r_Enew = (G * (M_sun + M_earth) * T_Enew**2 / (4 * np.pi**2))**(1/3)

    print(f"Desired Earth orbit in AU: {r_Enew.to(u.au):.3f}")

    # calculate the energy difference, i.e. the difference in potential energy
    E = G * M_sun * M_earth * ( 1/r_E - 1/r_Enew)

    print(f"Energy difference: {E.to(u.Joule):.1e}")

    # Caclute the power needed to provide this energy difference in 42 years
    P = (E / tau).to(u.Watt)

    print(f"Power needed: {P:.1e}")




    # INCLINATION CHANGE ------------------------------------------

    # Now what power would be needed to change the inclination of our orbit by 1 deg
    delta_i = 1 * u.deg

    # Earth's current speed in its orbit at apohelion, where the inclination 
    # change would be most efficient
    V = 29.29 * u.km / u.s 

    # delta V equation from Wikipedia, assuming eccentricity is zero
    delta_V = 2 * V * np.sin(delta_i / 2)

    print(f"Delta V needed: {delta_V:.1f}")

    # calculate the energy difference, i.e. the difference in kinetic energy
    E = 0.5 * M_earth * delta_V**2

    print(f"Energy difference: {E.to(u.Joule):.1e}")

    # Caclute the power needed to provide this energy difference in 42 years, that is 84 thrusts of 1h each
    P = (E / (84 * u.h)).to(u.Watt)

    print(f"Power needed for a {delta_i} inclination change: {P:.1e}")


