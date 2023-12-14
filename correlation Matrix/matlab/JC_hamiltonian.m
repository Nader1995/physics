
function [H_total] = JC_hamiltonian(B, BD, A, AD)
% Calculate the Hamiltonian 
% In each step, we call this function to calculate the total hamiltonian

% ======================================================================
% Constants
% According to the paper, this must hold: g * x_zpf << lambda,...
% Here I ignored the H_I and H_M, hence g = 0, Omega_M = 0

% lambda = 1 ;
% Omega = 1 ;          % For the photonic oscillation

% ======================================================================
% Defining Different Terms of Hamiltonian
% The Hamiltonian is time dependent

H_BS = (AD * B + BD * A) ;
H_C  = BD * B + AD * A ; 

H_total = H_BS + H_C ;

end
