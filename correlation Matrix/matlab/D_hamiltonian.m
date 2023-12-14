
function [H_total] = D_hamiltonian(C, B, BD, A, AD)
% Calculate the Hamiltonian 
% In each step, we call this function to calculate the total hamiltonian

% ======================================================================
% Constants

omega_A = 1 ;   % Oscillator A Energy
omega_B = 1 ;   % Oscillator B Energy
J = 1 ;         % Coupling Strength

% ======================================================================
% Defining Different Terms of Hamiltonian
% The Hamiltonian is time dependent

H_I = J * (AD * A * BD * B) ;

H_C  = omega_B * (C - 2 * BD * B) + omega_A * (AD * A) ; 

H_total = H_I + H_C ;

end

