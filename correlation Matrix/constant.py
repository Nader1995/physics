
# main calls matrixCalculator. It also receives the dimension of Fock state from user and updates the global constant N
# MatrixCalculator calculates W for our vectorised equation by calling dynamicCalculator,
# and comparing each time evolution by all the elements inside matrix C. So it also calls matrixGenerator.
# DynamicCalculator calls matrixGenerator and calculates the time evolution of each element.
# MatrixGenerator generates our main matrix C by receiving the dimension of fock space from user.
# Main passes the matrix W and constant N to matlab code.
