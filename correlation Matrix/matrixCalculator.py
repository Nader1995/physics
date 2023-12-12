import sympy as sp
import numpy as np
from matrixGenerator import cVec
from constant import N

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

preResult = []
for i in range(len(cVec)):

    eta = Commutator(HJC, cVec[i])
    preResult.append(normal_ordered_form(eta.doit().expand()))       # eta in the (separated) expanded form

print(preResult)

newList = []
for i in range(len(cVec)):

    newElement = 0
    for j in range(len(preResult[i].args)):

        newElement = newElement + preResult[i].args[j]
        for k in range(len(preResult[i].args[j].args)):

            if preResult[i].args[j].args[k] == a**2 or preResult[i].args[j].args[k] == Dagger(a)**2:

                newElement = newElement - preResult[i].args[j]
                break

    newList.append(newElement)

print(newList)
# print(etaSep.simplify())
# etaSep.args[1] = 0
# print(etaSep.args[1].args)
# print(len(etaSep.args[1].args))
# print(etaSep.args[1].args[5])

# result = sp.Mul(etaSep * a)
# print(etaSep)
# print(sp.srepr(etaSep))
# print(result.doit())
# print(normal_ordered_form(result.doit().expand()))
# complete form: pauli.qsimplify_pauli(normal_ordered_form(eta.doit().expand()))
# new = normal_ordered_form(eta.doit().expand())

# print(new.args[1])
# print(np.real(np.absolute(new.args[1])))
