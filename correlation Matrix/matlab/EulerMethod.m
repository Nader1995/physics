
% ******************************************
% *                                        *
% * August 22, 2023                        *
% * TED College                            *
% * Euler Method,                          *
% * Time evolution of spin population      *
% *                                        *
% * I find the time evolution of ...       *
% * ... sigma ^ dagger sigma using ...     *
% * ... using the famous Euler's method    *
% *                                        *
% *                                        *
% ******************************************\

tic ; 
clc ;
hold on

% =======================================================================
% Constants:

T = 2 * pi ;                % Totoal Time, (2 * pi) for a complete rotation 
step = 1/16 ;               % Powers of 2 : 16, 32, 64, ...
alpha = 5 ;                 % Degree of Coherency of alpha
dBoson = 4 ;               % Dimension of bosonic subspace
dSpin = 2 ;                 % Dimension of qubit subspace  

% =======================================================================
% Defining Annihilation and Creation Operators 

bI = eye(dBoson) ;
bVec = sqrt(1:dBoson - 1) ;
b = diag(bVec, 1) ;
bd = b' ;

sI = eye(dSpin) ;
sVec = sqrt(1:dSpin - 1) ;
s = diag(sVec, 1) ;
sd = s' ;

% =======================================================================
% Defining Annihilation and Creation Operators of Bosonic and Atomic Fields
% The Order of Subspaces is as follows:
% Atomic + Bosonic = S + B

S = kron(s, bI) ;
B = kron(sI, b) ;

SD = S';
BD = B';

% =======================================================================
% Defining the States

bket = zeros(dBoson, 1) ;

bketZero = bket ;
bketZero(1) = 1;        % |0 > : ground state, [1, 0, 0, ...]

bketOne = bket ;
bketOne(2) = 1 ;        % |1 > : first excited state: [0, 1, 0, ...]

bketTwo = bket ;
bketTwo(3) = 1 ;        % |2 > : ground state, [0, 0, 1, ...]

bketThree = bket ;
bketThree(4) = 1 ;      % |3 > : ground state, [0, 0, 0, 1, ...]

sket = zeros(dSpin, 1) ;

sketZero = sket ;
sketZero(1) = 1 ;       % |0 > : ground state, [1, 0, 0, ...]

sketOne = sket ;
sketOne(2) = 1 ;        % |1 > : first excited state: [0, 1, 0, ...]

ketAlpha = exp(-0.5 * abs(alpha) ^ 2) * expm(alpha * bd) * expm(-conj(alpha) * b) * bketZero ;     % coherent state alpha

% state = kron(sketZero * sketZero', ketAlpha * ketAlpha') ;
state = kron(sketZero * sketZero', bketOne * bketOne'); 

% =======================================================================
% main part:

Time = 0 : step : T ;

hop = zeros((1/step) * (round(T) + 1), 1) ;
hopd = zeros((1/step) * (round(T) + 1), 1) ;
pop = zeros((1/step) * (round(T) + 1), 1) ;
sigma = zeros((1/step) * (round(T) + 1), 1) ;

counter = 1 ;

hop(counter) = trace(state * S * BD) ;
hopd(counter) = trace(state * SD * B) ;
pop(counter) = trace(state * BD * B) ;
sigma(counter) = trace(state * SD * S) ;

for time = 0 : step : T 
    
    sigma(counter + 1) = sigma(counter) + 1i * step * (hop(counter) - hopd(counter)) ;
    
    pop(counter + 1) = pop(counter) + 1i * step * (-hop(counter) + hopd(counter)) ;
    
    hop(counter + 1) = hop(counter)  + 1i * step * (sigma(counter) + 2 * pop(counter) * sigma(counter) - pop(counter)) ;
    
    hopd(counter + 1) = hopd(counter) - 1i * step * (sigma(counter) + 2 * pop(counter) * sigma(counter) - pop(counter)) ;
    
    disp(sigma(counter)) ;
    
    counter = counter + 1;
    
end

plot(sigma);
