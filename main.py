# ********************************************************************
# *                                                                  *
# * correlation Matrix, v 1.00                                       *
# * Date: 07.Dec.2023                                                *
# * Nader Safari,                                                    *
# * Copyleft 🄯 2023 N. Safari. All left reserved!                    *
# *                                                                  *
# ********************************************************************

# Source: https://qiskit.org/textbook/ch-quantum-hardware/cQED-JC-SW.html

import numpy as np

from constant import N
from matrixGenerator import cVec
from matrixCalculator import W, wTime, polishTime, commuteTime

print(N)
print(cVec)
print(W)

print('commuteTime: ', commuteTime)
print('PolishTime: ', polishTime)
print('wTime: ', wTime)

np.savetxt('W.mat', W)
np.savetxt('N.mat', [N])
