
% *********************************************** 
% *                                             *
% * 15 November, 2023                           *
% * TED College                                 *
% * Classical Entanglement Limit, part I        *
% *                                             *
% * |n> = 1                                     *
% *                                             *
% * Here I study the dynamics                   *
% * of qubit-cavity system                      *
% * For "Jaynes-Cummings Hamiltonian"           *
% *                                             *
% * I used the classical entanglement method,   *
% * one can cancel some higher order terms      *
% * in the matrix of moments, considering some  *
% * specific initial states                     *
% *                                             *
% ***********************************************

tic 
clc ;   
hold on ;

% =======================================================================
% Constants:

T = 2 * pi ;            % Totoal Time, (2 * pi) for a complete rotation 
step = 1/32 ;           % Powers of 2 : 16, 32, 64, ...

dimension = 3 ;        % Dimension of bosonic subspace
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

C = cell(16, 1);

C{1} = kron(qI, aI) ;
C{2} = kron(qI, a) ;
C{3} = kron(sigma, aI) ;
C{4} = kron(sigma, a) ;
C{5} = kron(qI, ad) ;
C{6} = kron(qI, ad * a) ;
C{7} = kron(sigma, ad) ;
C{8} = kron(sigma, ad * a) ;
C{9} = kron(sigmad, aI) ;
C{10} = kron(sigmad, a) ;
C{11} = kron(sigmad * sigma, aI) ;
C{12} = kron(sigmad * sigma, a) ;
C{13} = kron(sigmad, ad) ;
C{14} = kron(sigmad, ad * a) ;
C{15} = kron(sigmad * sigma, ad) ;
C{16} = kron(sigmad * sigma, ad * a) ;

% =======================================================================
% Defining the States

% First of all, let's define our pure states 
ket = zeros(dimension,1) ;
ketZero = ket ;
ketZero(1) = 1;          % |0 > : ground state, [1, 0, 0, ...]
ketOne = ket ;
ketOne(2) = 1 ;          % |1 > : first excited state: [0, 1, 0, ...]
ketTwo = ket ;
ketTwo(3) = 1 ;          % |2 > : ground state, [0, 0, 1, ...]

qKet = zeros(qDimension, 1) ;
qKetZero = qKet ;
qKetZero(1) = 1 ;       % |0 > : ground state, [1, 0, 0, ...]
qKetOne = qKet ;
qKetOne (2) = 1 ;

state = kron(qKetZero * qKetZero', ketOne * ketOne'); 
% state = kron(qKetOne * qKetOne', ketZero * ketZero'); 
% state = kron(qKetZero * qKetZero', ketTwo * ketTwo'); 

C_ini = zeros(16, 1) ;

for i = 1:16
        
    C_ini(i) = trace(state * C{i}) ;
end

% =======================================================================
% Lindbladian
X = zeros(16, 16) ;

X(2, 3) = -1 ;

X(3, 2) = -1 ;
X(3, 12) = 2 ;

X(5, 9) = +1 ;

X(6, 7) = -1 ;
X(6, 10) = +1 ;

X(7, 11) = +1 ;
X(7, 16) = +2 ;
X(7, 6) = -1 ;

X(8, 12) = +1 ;

X(9, 5) = +1 ;
X(9, 15) = -2 ;

X(10, 11) = -1 ;
X(10, 16) = -2 ;
X(10, 6) = +1 ;

X(11, 7) = +1 ;
X(11, 10) = -1 ;

X(12, 8) = +1 ;

X(14, 15) = -1 ;

X(15, 14) = -1 ;

% To calculate C2
lindbladian = 1i * X ;

% ======================================================================
% Main: Integral Method

Time = 0 : step : T ;   % total steps, we need this in our plot

resultC = zeros(round((1/step) * T), 1) ;     % C2 is stored here
  
counter = 1 ;
  
for time = 0:step:T
   
   co = expm(lindbladian * time) ;
   result = co * C_ini ;
  
   % Calculations regarding C2
%    resultC(counter, 1) = result(11, 1) ;        % for sigmad * sigma
   resultC(counter, 1) = result(6, 1) ;        % for ad * a
   
   disp( " ================================================ " );
   disp( resultC(counter, 1) );
 
%    disp("I am here: ");
%    disp(time);
    
   counter = counter + 1;
    
end

% =======================================================================
% Plot the Result                         

plot(Time, real(resultC));

elapsed_Time = toc ;
disp('Run Iime (minute):') ;
disp(elapsed_Time / 60);
