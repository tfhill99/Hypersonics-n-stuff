# Hypersonics-n-stuff

Works to calculate the entry trajectory of a re-entry capsule given initial conditions.
Then calculates the acceleration, pressure, lift and drag coefficients for a specified geometry. 

Finally, iteratively converges to a C_d, C_l and non-optimized flight path trajectory.

Reads in an optimized database to calculate q_w using Tauber Menees at stagnation and along the body.

Calculates the 1D heat transfer solution for various material thicknesses and optimizes a TPS for a payload back temperature of 70C.
