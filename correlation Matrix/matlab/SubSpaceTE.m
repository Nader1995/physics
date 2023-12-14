
% *********************************************** 
% *                                             *
% * 3 September, 2023                           *
% * TED College                                 *
% * Sub-space Time Evolution,                   *
% *                                             *
% *                                             *
% * Here I study the dynamics                   *
% * of qubit-cavity system                      *
% * For "dispersive Hamiltonian"                *
% *                                             *
% * I used the subspace-time derivation method  *
% * also known as SubSpace TE method            *
% * Here, matrix C is based on                  *
% * My own approach                             *
% *                                             *
% * Dispersive H                                *
% ***********************************************

tic 
clc ;   
hold on ;

% =======================================================================
% Constants:

T = 2 * pi ;            % Totoal Time, (2 * pi) for a complete rotation 
step = 1/32 ;           % Powers of 2 : 16, 32, 64, ...
alpha = 5 ;             % Degree of Coherency of alpha

% dimension of bosonic field is approved to be the best fit for a ...
% ... coherent state with value of alpha = 5
dimension = 68 ;        % Dimension of bosonic subspace
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
% Defining Annihilation and Creation Operators of Bosonic and Atomic Fields
% The Order of Subspaces is as follows:
% Atomic + Bosonic = Q + A 

I = kron(qI, aI) ;
Q = kron(sigma, aI) ;
A = kron(qI, a) ;

QD = Q';            % QD = sigma^dagger
AD = A';            % A = a^dagger

ground = Q * QD ;   % Atomic ground state
excited = QD * Q ;  % Atomic excited state

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
ketThr = ket ;
ketThree(4) = 1 ;        % |3 > : ground state, [0, 0, 0, 1, ...]

qKet = zeros(qDimension, 1) ;
qKetZero = qKet ;
qKetZero(1) = 1 ;       % |0 > : ground state, [1, 0, 0, ...]
qKetOne = qKet ;
qKetOne (2) = 1 ;

ketAlpha = exp(-0.5 * abs(alpha) ^ 2) * expm(alpha * ad) * expm(-conj(alpha) * a) * ketZero ;     % coherent state alpha

state = kron(qKetZero * qKetZero', ketAlpha * ketAlpha') ;
% state = kron(qKetZero * qKetZero', ketOne * ketOne'); 
% state = kron(qKetZero * qKetZero', ketTwo * ketTwo'); 
% state = kron(qKetOne * qKetOne', ketZero * ketZero'); 
% state = kron(qKetOne * qKetOne', ketOne * ketOne'); 

atom = cell(2, 2);
atom{1, 1} = ground ;
atom{1, 2} = Q ;
atom{2, 1} = QD ;
atom{2, 2} = excited ;


field = cell(2, 2);
field{1, 1} = I ;
field{1, 2} = A ;
field{2, 1} = AD ;
field{2, 2} = AD * A ;

C = cell(16, 1);

C{1} = atom{1} * field{1} ;
C{2} = atom{1} * field{2} ;
C{3} = atom{2} * field{1} ;
C{4} = atom{2} * field{2} ;
C{5} = atom{1} * field{3} ;
C{6} = atom{1} * field{4} ;
C{7} = atom{2} * field{3} ;
C{8} = atom{2} * field{4} ;
C{9} = atom{3} * field{1} ;
C{10} = atom{3} * field{2} ;
C{11} = atom{4} * field{1} ;
C{12} = atom{4} * field{2} ;
C{13} = atom{3} * field{3} ;
C{14} = atom{3} * field{4} ;
C{15} = atom{4} * field{3} ;
C{16} = atom{4} * field{4} ;

C_ini = zeros(16, 1) ;

for i = 1:16
        
    C_ini(i) = trace(state * C{i}) ;
end

% =======================================================================
% Lindbladian
X = zeros(16, 16) ;

X(2, 2) = -1 ;
X(4, 4) = -2 ;
X(5, 5) = +1 ;
X(6, 1) = +1 ;
X(7, 7) = +2 ;
X(8, 3) = +2 ;
X(10, 10) = -1 ;
X(12, 12) = -2 ;
X(13, 13) = +1 ;
X(14, 9) = +1 ;
X(15, 15) = +2 ;
X(16, 11) = +2 ;

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
   resultC(counter, 1) = result(6, 1) + result(16, 1) ;
   
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


