import sympy as sp
import numpy as np

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

sd = pauli.SigmaMinus()      # sd = sPlus = Dagger(sigma)
s = pauli.SigmaPlus()      # s = sMinus = sigma

# Define J-C Hamiltonian
HJC = sp.I * (ad * s + a * sd)       # sympy.I = √ (-1) = i

eta = Commutator(HJC, a * s)

# complete form: pauli.qsimplify_pauli(normal_ordered_form(eta.doit().expand()))
new = normal_ordered_form(eta.doit().expand())

print(new)
# print(new.args[1].args[1])
# print(normal_ordered_form(new.args[1] / (sp.I * a * s * sd)))
# print(abs(new.args[1]) == abs(sp.I * a * s * sd))
# print(abs(new.args[1]).doit() == abs(sp.I * a * s * sd))

# print(np.real(np.absolute(new.args[1])))

term = 0
term = new.args[1]
anotherTerm = new.args[1].as_coefficients_dict()
# print(term * anotherTerm.values())
print(new.args[1])
print(anotherTerm.values())

for i in anotherTerm.keys():
    print(i)
