

% *********************************************** 
% *                                             *
% * 16 December, 2023                           *
% * TED College                                 *
% * Python-MatLab Cooperation                   *
% *                                             *
% * |n> = n                                     *
% *                                             *
% * Here I study the dynamics                   *
% * of qubit-cavity system                      *
% * For "Jaynes-Cummings Hamiltonian"           *
% *                                             *
% * We find W usin the python code, and         *
% * import it here                              *
% *                                             *
% ***********************************************

tic 
clc ;   
hold on ;

% To import matrix W and integer N from python code
W = load('W.mat', '-ASCII') ;
N = load('N.mat', '-ASCII') ;

% =======================================================================
% Constants:

T = 10 * pi ;            % Totoal Time, (2 * pi) for a complete rotation 
step = 1/32 ;           % Powers of 2 : 16, 32, 64, ...
alpha = 4 ;

dimension = N ;        % Dimension of bosonic subspace
qDimension = 2 ;        % Dimension of qubit subspace  

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
% Defining the States

% First of all, let's define our pure states 
ket = zeros(dimension,1) ;
ketZero = ket ;
ketZero(1) = 1 ;          % |0 > : ground state, [1, 0, 0, ...]
ketOne = ket ;
ketOne(2) = 1 ;          % |1 > : first excited state: [0, 1, 0, ...]
ketTwo = ket ;
ketTwo(3) = 1 ;          % |2 > : ground state, [0, 0, 1, ...]
ketThree = ket ;
ketThree(4) = 1 ;          % |3 > : ground state, [0, 0, 0, 1, ...]
ketFour = ket ;
ketFour(5) = 1 ;          % |4 > : ground state, [0, 0, 0, 0, 1, ...]
ketFive = ket ;
ketFive(6) = 1 ;          % |5 > 
ketSix = ket ;
ketSix(7) = 1 ;          % |6 > 
ketSeven = ket ;
ketSeven(8) = 1 ;          % |7 > 
ketEight = ket ;
ketEight(9) = 1 ;          % |8 >
ketNine = ket ;
ketNine(10) = 1 ;          % |9 >
ketTen = ket ;
ketTen(11) = 1 ;          % |10 >

ketAlpha = exp(-0.5 * abs(alpha) ^ 2) * expm(alpha * ad) * expm(-conj(alpha) * a) * ketZero ;     % coherent state alpha

qKet = zeros(qDimension, 1) ;
qKetZero = qKet ;
qKetZero(1) = 1 ;       % |0 > : ground state, [1, 0, 0, ...]
qKetOne = qKet ;
qKetOne (2) = 1 ;


stateSeven = ketSeven * ketSeven' ;
stateSix = ketSix * ketSix' ;
stateFive = ketFive * ketFive' ;
stateFour = ketFour * ketFour' ;
stateThree = ketThree * ketThree' ;
stateTwo = ketTwo * ketTwo' ;
stateOne = ketOne * ketOne' ;
coherentState = ketAlpha * ketAlpha' ;

qState = qKetZero * qKetZero' ;
state = coherentState;
% =======================================================================
% Initial values

S = zeros(2, 2) ;

S(1, 1) = trace(qState * sigma * sigmad) ;
S(1, 2) = trace(qState * sigma) ;
S(2, 1) = trace(qState * sigmad) ;
S(2, 2) = trace(qState * sigmad * sigma) ;

Ad = cell(N, 1) ;

for i=1:N
    
    Ad{i} = ad^(i-1) ;
end

A = cell(1, N) ;

for i=1:N
    
    A{i} = a^(i-1) ;
end

AdA = cell(N, N) ;

for i=1:N
    for j=1:N
        
        AdA{i,j} = Ad{i} * A{j} ;
    end
end

C = zeros(N, N) ;

for i=1:N
    for j=1:N
        
%         C(i, j) = trace(state * AdA{i,j}) ;
%         When we are having coherent State
        C(i, j) = trace(state * AdA{i,j})/(alpha^2) ;
        
    end
end

% ==================================================
% Vectorization made easy

coMatrix = kron(S, C) ;
coMatrixTrans = coMatrix' ;
coVector = coMatrixTrans(:) ;

% ======================================================================
% Main: Integral Method

Time = 0 : step : T ;   % total steps, we need this in our plot

resultC = zeros(round((1/step) * T), 1) ;     % C2 is stored here
  
counter = 1 ;
  
for time = 0:step:T
   
   co = expm(1i * W * time) ;
   result = co * coVector ;
   
   resultC(counter, 1) = result(2*N + 2, 1) + result(2 * N^2 + 3*N + 2, 1) ; 
   
   disp( " ================================================ " );
   disp( resultC(counter, 1) );
 
   disp("I am here: ");
   disp(counter);
    
   counter = counter + 1;
    
end

% =======================================================================
% Plot the Result                         

plot(Time, real(resultC));

elapsed_Time = toc ;
disp('Run Iime (minute):') ;
disp(elapsed_Time / 60);




