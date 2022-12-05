# modules
import numpy as np

# Constants and Parameters
# ========================
kB = float(1.38E-23)            # Boltzmann constant [J/K]
eV = float(1.6E-19)             # ElectroVolt [J]
Da = float(1.6E-27)             # The mass of an atom (atomic mass) is often expressed in the Dalton unit.
# According to the definition, 1 Dalton is defined as ¹/₁₂ of the mass of
# a single ¹²C atom at rest [kg].

mAr = float(40 * Da)            # Mass of Argon atom [kg]
sigma = float(3.40e-10)         # σ of the Lennard-Jones potential; σ = 3.40 Å [m].
epsilon = float(117 * kB)       # ɛ of the Lennard-Jones potential [J].

# <===  Units in the reduced units system  ===>
ul = sigma                      # unit of length; l*
um = mAr                        # unit of length; m*
ut = np.sqrt(mAr*np.sqrt(sigma)/epsilon)       # unit of time; t*
uE = epsilon                                   # unit of energy; ɛ
uv = np.sqrt(epsilon/mAr)                      # unit of velocity; v*
uT = epsilon/kB                                # unit of temperature; T*

# <===  MD parameters  ===>
N = int(100)                     # Number of particles (Argon atoms)
L = int(30)                     # The length of the simulation box [l*]
h = float(1/128)                # It is better to choose a power of 2 for the time step [t*].
T0 = float(10)                  # Initial temperature [T*]
t_max = float(10)               # The maximum time of simulation
# ------------------------------------------
