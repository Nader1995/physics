import sympy as sp
from constant import N

# import operator relations and define them
from sympy import Matrix
from sympy.physics.quantum.boson import BosonOp
from sympy.physics.quantum import pauli, Dagger, TensorProduct
from sympy.physics.quantum.operatorordering import normal_ordered_form

# It is supposed to make outcome pretty
sp.init_printing(use_unicode=False, wrap_line=False)

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

sd = pauli.SigmaMinus()  # sd = sPlus = Dagger(sigma)
s = pauli.SigmaPlus()  # s = sMinus = sigma

Sigma = Matrix([[s * sd, s],
                [sd, sd * s]])


def f(x, y):
    return a ** y


A = Matrix(1, N, f)


def g(x, y):
    return ad ** x


Ad = Matrix(N, 1, g)

Alpha = Ad * A

C = normal_ordered_form(TensorProduct(Sigma, Alpha))

cVec = []
k = 0
for i in range(2*N):
    for j in range(2*N):
        cVec.append(C[i, j])
        k += 1
