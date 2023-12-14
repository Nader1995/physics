import sympy as sp
from matrixGenerator import cVec
from constant import N

# import operator relations and define them
from sympy.matrices import zeros
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
HJC = ad * s + a * sd       # sympy.I = √ (-1) = i

preResult = []
for i in range(len(cVec)):

    eta = Commutator(HJC, cVec[i])
    preResult.append(normal_ordered_form(eta.doit().expand()))       # eta in the (separated) expanded form

# ==============================================================
# Stage I: cleaning the correlation matrix of all unwanted terms

newList = []
for i in range(len(cVec)):

    newElement = 0
    for j in range(len(preResult[i].args)):

        newElement = newElement + preResult[i].args[j]
        createdElement = 1
        for k in range(len(preResult[i].args[j].args)):

            createdElement = createdElement * preResult[i].args[j].args[k]
            if preResult[i].args[j].args[k] == a**N or preResult[i].args[j].args[k] == Dagger(a)**N:

                newElement = newElement - preResult[i].args[j]
                break

            # These two if-clauses are Hamiltonian dependent: REMEMBER
            # HJC = sp.I * (ad * sd * s * sd + a * s * sd * s) & Commutator(HJC, a * s)
            if preResult[i].args[j].args[k] == sd and len(preResult[i].args[j].args)-k == 3:

                newElement = newElement - preResult[i].args[j] + createdElement
                break

            elif preResult[i].args[j].args[k] == s and len(preResult[i].args[j].args)-k == 3:

                newElement = newElement - preResult[i].args[j] + createdElement
                break

    newList.append(newElement)

result = newList

# print(preResult)
# print(result)

# ==================================================
# Stage II: Calculating W

W = zeros(N**2, N**2)

for i in range(len(result)):

    for j in range(len(result[i].args)):

        print(result[i].args[j].args)
        for k in range(len(cVec)):

            if abs(result[i].args[j]) == abs(cVec[k]):

                print(True)
