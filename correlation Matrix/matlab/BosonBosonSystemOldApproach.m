
% **********************************
% *                                *
% * July 11, 2023                  *
% * TED College                    *
% * Boson Boson Model,             *
% * The old approach 1st Version   *
% *                                *
% * Here I study the dynamics      *
% * of a two-level Boson system,   *
% * interacting with another       *
% * bosonic system.                *
% *                                *
% *                                *
% **********************************

tic 
clc ;   
hold on ;

% =======================================================================
% Constants:

T = 2 * pi ;          % Totoal Time, (2 * pi) for a complete rotation 
step = 1/16 ;        % Powers of 2 : 16, 32, 64, ...
alpha = 5 ;        % Degree of Coherency of alpha
dimension1 = 60 ;       % Dimension of bosonic subspace
dimension2 = 2 ;       % Dimension of qubit subspace  

% =======================================================================
% Defining Annihilation and Creation Operators 

aI = eye(dimension1) ;
aVec = sqrt(1:dimension1 - 1) ;
a = diag(aVec, 1) ;
ad = a' ;

bI = eye(dimension2) ;
bVec = sqrt(1:dimension2 - 1) ;
b = diag(bVec, 1) ;
bd = b' ;

% =======================================================================
% Defining Annihilation and Creation Operators of Bosonic and Atomic Fields
% The Order of Subspaces is as follows:
% Bosonic + Atomic = A + B

A = kron(a, bI) ;
B = kron(aI, b) ;

AD = A';
BD = B';

% =======================================================================
% Defining the Hamiltonian for initial time

H_total = JC_hamiltonian(A, AD, B, BD) ;

% =======================================================================
% Defining the States

% First of all, let's define our pure states 
aket = zeros(dimension1,1) ;
aketZero = aket ;
aketZero(1) = 1;        % |0 > : ground state, [1, 0, 0, ...]
aketOne = aket ;
aketOne(2) = 1 ;        % |1 > : first excited state: [0, 1, 0, ...]
aketTwo = aket ;
aketTwo(3) = 1 ;         % |2 > : ground state, [0, 0, 1, ...]
aketThree = aket ;
aketThree(4) = 1 ;       % |3 > : ground state, [0, 0, 0, 1, ...]

bket = zeros(dimension2, 1) ;
bketZero = bket ;
bketZero(1) = 1 ;       % |0 > : ground state, [1, 0, 0, ...]

ketAlpha = exp(-0.5 * abs(alpha) ^ 2) * expm(alpha * ad) * expm(-conj(alpha) * a) * aketZero ;     % coherent state alpha

state = kron(ketAlpha * ketAlpha', bketZero * bketZero') ;
% state = kron(aketTwo * aketTwo', bketZero * bketZero'); 

% =======================================================================
% main:

Time = 0 : step : T ;
result1 = zeros((1/step) * (round(T) + 1), 1) ;
j = 1 ;

for t = 0:step:T
    U = expm(-1i * H_total * t);
    
    N = U' * (BD * B) * U ;

    result1(j, 1) = trace(state * N);
  
    disp("I'm Here: ") ;
    disp(j) ;
    disp(result1(j, 1));
    
    j = j + 1 ;    
end

while result1(end) == 0
    result1(end) = [] ;
end

% =======================================================================
% Plot the Result
% final_Result1 is for HH case: Indistinguishable

final_Result1 = real(transpose(result1)) ;

plot(Time, final_Result1) ; % HH
elapsed_Time = toc ;
disp('Run Iime (minute): ') ;
disp(elapsed_Time / 60);
