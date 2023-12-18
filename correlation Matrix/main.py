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
from matrixCalculator import W, wTime, polishTime, commuteTime, result, preResult

print(N)
print(cVec)
print(preResult)
print(result)

# The time it takes to calculate the commutation relation (This is the longest)
print('commuteTime: ', commuteTime)
# The time it takes to delete all unnecessary elements
print('PolishTime: ', polishTime)
# The time it takes to create W matrix
print('wTime: ', wTime)

# To save as text and load it from MatLab
np.savetxt('W.mat', W)
# Since it does not save integers, save it as 1-D array [N]
np.savetxt('N.mat', [N])
