
% *************************************************
% *                                               *
% * June 28, 2023                                 *
% * TED College                                   *
% * Spin Boson Model,                             *
% * The old approach 1st Version                  *
% *                                               *
% * Here I study the dynamics                     *
% * of Atomic-Bosonic system                      *
% *                                               *
% * Here, I used the analytical...                *
% * ... approach (Solving Hisenberg Equation)     *
% *                                               *
% *************************************************

tic 
clc ;   
hold on ;

% =======================================================================
% Constants:

T = 2 * pi ;          % Totoal Time, (2 * pi) for a complete rotation 
step = 1/32 ;        % Powers of 2 : 16, 32, 64, ...
alpha = 5 ;        % Degree of Coherency of alpha
dimension = 4;       % Dimension of bosonic subspace: 68 for best fit
qDimension = 2 ;       % Dimension of qubit subspace  

% =======================================================================
% Defining Annihilation and Creation Operators 

aI = eye(dimension) ;
vec = sqrt(1:dimension - 1) ;
a = diag(vec, 1) ;
ad = a' ;

qI = eye(qDimension) ;
qVec = sqrt(1:qDimension - 1) ;
sigma = diag(qVec, 1) ;
sigmad = sigma' ;

% =======================================================================
% Defining Annihilation and Creation Operators of Bosonic and Atomic Fields
% The Order of Subspaces is as follows:
% Atomic + Bosonic = Q + A 

I = kron(qI, aI) ;
Q = kron(sigma, aI) ;
A = kron(qI, a) ;

QD = Q';
AD = A';

% =======================================================================
% Defining the Hamiltonian for initial time

H_total = JC_hamiltonian(Q, QD, A, AD) ;
% H_total = D_hamiltonian(I, Q, QD, A, AD) ;

% =======================================================================
% Defining the States

% First of all, let's define our pure states 
ket = zeros(dimension,1) ;
ketZero = ket ;
ketZero(1) = 1;         % |0 > : ground state, [1, 0, 0, ...]
ketOne = ket ;
ketOne(2) = 1 ;         % |1 > : first excited state: [0, 1, 0, ...]
ketTwo = ket ;
ketTwo(3) = 1 ;         % |2 > : ground state, [0, 0, 1, ...]
ketThree = ket ;
ketThree(4) = 1 ;       % |3 > : ground state, [0, 0, 0, 1, ...]

qKet = zeros(qDimension, 1) ;
qKetZero = qKet ;
qKetZero(1) = 1 ;       % |0 > : ground state, [1, 0, 0, ...]
qKetOne = qKet ;
qKetOne (2) = 1 ;

ketAlpha = exp(-0.5 * abs(alpha) ^ 2) * expm(alpha * ad) * expm(-conj(alpha) * a) * ketZero ;     % coherent state alpha

% state = kron(qKetZero * qKetZero', ketAlpha * ketAlpha') ;
% state = kron(qKetZero * qKetZero', ketOne * ketOne'); 
state = kron(qKetZero * qKetZero', ketTwo * ketTwo');
% state = kron(qKetZero * qKetZero', ketThree * ketThree'); 
% state = kron(qKetOne * qKetOne', ketZero * ketZero'); 
% state = kron(qKetOne * qKetOne', ketOne * ketOne'); 

% =======================================================================
% main:

Time = 0 : step : T ;
result1 = zeros(round((1/step) * T), 1);     % calculated result is stored here 
j = 1 ;

for t = 0:step:T
    U = expm(-1i * H_total * t);
    
%     N = U' * (QD * Q) * U ;
    N = U' * (AD * A) * U ; 

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

