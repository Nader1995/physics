
% *********************************************** 
% *                                             *
% * 6 September, 2023                           *
% * TED College                                 *
% * Sub-space Time Evolution,                   *
% *                                             *
% *                                             *
% * Here I study the dynamics                   *
% * of qubit-qubit system                       *
% * For "Jaynes-Cummings Hamiltonian"           *
% *                                             *
% * I used the subspace-time derivation method  *
% * also known as SubSpace TE method            *
% *                                             *
% *                                             *
% ***********************************************

tic 
clc ;   
hold on ;

% =======================================================================
% Constants:

T = 2 * pi ;            % Totoal Time, (2 * pi) for a complete rotation 
step = 1/32 ;           % Powers of 2 : 16, 32, 64, ...

% We don't get into coherent state example becuase the dimensions are low
% alpha = 5 ;             % Degree of Coherency of alpha

dimension = 2 ;        % Dimension of bosonic subspace

% =======================================================================
% Defining Annihilation and Creation Operators 
% The Order of Subspaces is as follows:
% Atom One + Atom Two = Sigma1 + Sigma2
% We are interested in population of the first system

I = eye(dimension) ;
vec = sqrt(1:dimension - 1) ;

sigma1 = diag(vec, 1) ;
sigma2 = diag(vec, 1) ;

Sigma1 = kron(sigma1, I) ;
Sigma2 = kron(I, sigma2) ;

% =======================================================================
% Defining the Hamiltonian for initial time

% H_total = JC_hamiltonian(Sigma1, Sigma1', Sigma2, Sigma2') ;
H_total = D_hamiltonian(kron(I, I), Sigma1, Sigma1', Sigma2, Sigma2') ;

% =======================================================================
% Defining the States

% First of all, let's define our pure states 
ket = zeros(dimension,1) ;
ketZero = ket ;
ketZero(1) = 1;         % |0 > : ground state, [1, 0, 0, ...]
ketOne = ket ;
ketOne(2) = 1 ;         % |1 > : first excited state: [0, 1, 0, ...]

% state = kron(ketZero * ketZero', ketOne * ketOne'); 
state = kron(ketOne * ketOne', ketZero * ketZero'); 
% state = kron(ketOne * ketOne', ketOne * ketOne'); 

% =======================================================================
% main:

Time = 0 : step : T ;
result1 = zeros(round((1/step) * T), 1);     % calculated result is stored here 
j = 1 ;

for t = 0:step:T

    U = expm(-1i * H_total * t);
    
    N = U' * (Sigma1' * Sigma1) * U * (Sigma2' * Sigma2) + U' * (Sigma1' * Sigma1) * U * (Sigma2 * Sigma2') ;
%     N = U' * (Sigma1' * Sigma1) * U * (Sigma2' * Sigma2) ;
%     N = U' * (Sigma1' * Sigma1) * U * (Sigma2 * Sigma2') ;
    
    result1(j, 1) = trace(state * N);
  
    disp("I'm Here: ") ;
    disp(j) ;
    disp(result1(j, 1));
    
    j = j + 1 ;    
end

% =======================================================================
% Plot the Result

final_Result1 = real(transpose(result1)) ;

plot(Time, final_Result1) ; 
elapsed_Time = toc ;
disp('Run Iime (minute): ') ;
disp(elapsed_Time / 60)

