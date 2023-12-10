# ********************************************************************
# *                                                                  *
# * correlation Matrix, v 1.00                                       *
# * Date: 07.Dec.2023                                                *
# * Nader Safari,                                                    *
# * Copyleft 🄯 2023 N. Safari. All left reserved!                    *
# *                                                                  *
# ********************************************************************

# Source: https://qiskit.org/textbook/ch-quantum-hardware/cQED-JC-SW.html

import sympy as sp
import numpy as np
from matrixGenerator import *
from matrixCalculator import *

# import operator relations and define them
from sympy.physics.quantum.boson import BosonOp
from sympy.physics.quantum import pauli, Dagger, Commutator
from sympy.physics.quantum.operatorordering import normal_ordered_form

# It is supposed to make outcome pretty
sp.init_printing(use_unicode=True)

# Annihilation operator
a = BosonOp('a')
ad = Dagger(a)

# Pauli matrices
sx = pauli.SigmaX()
sy = pauli.SigmaY()
sz = pauli.SigmaZ()

# Qubit raising and lowering operators,
# NOTICE the inverse relation: sz = I - 2 * pauli.SigmaMinus() * pauli.SigmaPlus(),
# sd = sPlus = pauli.SigmaMinus() = Dagger(sigma)
# s = sMinus = pauli.SigmaPlus() = sigma

sd = pauli.SigmaMinus()      # sPlus = Dagger(sigma) = sd
s = pauli.SigmaPlus()      # sMinus = sigma = s

# Define J-C Hamiltonian
HJC = sp.I * (ad * s + a * sd)       # sympy.I = √ (-1) = i

eta = Commutator(HJC, s)

# complete form: pauli.qsimplify_pauli(normal_ordered_form(eta.doit().expand()))
new = normal_ordered_form(eta.doit().expand())

print(new)
print(np.real(np.absolute(new.args[1])))
