
% *********************************************** 
% *                                             *
% * 31 November, 2023                           *
% * TED College                                 *
% * Classical Entanglement Limit, part II       *
% *                                             *
% * |n> = 2                                     *
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

C = cell(36, 1);

C{1} = kron(qI, aI) ;
C{2} = kron(qI, a) ;
C{3} = kron(qI, a * a) ;
C{4} = kron(sigma, aI) ;
C{5} = kron(sigma, a) ;
C{6} = kron(sigma, a * a) ;
% =================================

C{7} = kron(qI, ad) ;
C{8} = kron(qI, ad * a) ;
C{9} = kron(qI, ad * a * a) ;
C{10} = kron(sigma, ad) ;
C{11} = kron(sigma, ad * a) ;
C{12} = kron(sigma, ad * a * a) ;
% =================================

C{13} = kron(qI, ad * ad) ;
C{14} = kron(qI, ad * ad * a) ;
C{15} = kron(qI, ad * ad * a * a) ;
C{16} = kron(sigma, ad * ad) ;
C{17} = kron(sigma, ad * ad * a) ;
C{18} = kron(sigma, ad * ad * a * a) ;
% =================================

C{19} = kron(sigmad, aI) ;
C{20} = kron(sigmad, a) ;
C{21} = kron(sigmad, a * a) ;
C{22} = kron(sigmad * sigma, aI) ;
C{23} = kron(sigmad * sigma, a) ;
C{24} = kron(sigmad * sigma, a * a) ;
% =================================

C{25} = kron(sigmad, ad) ;
C{26} = kron(sigmad, ad * a) ;
C{27} = kron(sigmad, ad * a * a) ;
C{28} = kron(sigmad * sigma, ad) ;
C{29} = kron(sigmad * sigma, ad * a) ;
C{30} = kron(sigmad * sigma, ad * a * a) ;
% =================================

C{31} = kron(sigmad, ad * ad) ;
C{32} = kron(sigmad, ad * ad * a) ;
C{33} = kron(sigmad, ad * ad * a * a) ;
C{34} = kron(sigmad * sigma, ad * ad) ;
C{35} = kron(sigmad * sigma, ad * ad * a) ;
C{36} = kron(sigmad * sigma, ad * ad * a * a) ;
% =================================

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

alpha = 0.95 ;
ketAlpha = exp(-0.5 * abs(alpha) ^ 2) * expm(alpha * ad) * expm(-conj(alpha) * a) * ketZero ;     % coherent state alpha

% state = kron(qKetZero * qKetZero', ketAlpha * ketAlpha') ;
% state = kron(qKetZero * qKetZero', ketOne * ketOne'); 
% state = kron(qKetOne * qKetOne', ketZero * ketZero'); 
state = kron(qKetZero * qKetZero', ketTwo * ketTwo'); 

C_ini = zeros(36, 1) ;

for i = 1:36
        
    C_ini(i) = trace(state * C{i}) ;
end

% =======================================================================
% Lindbladian
X = zeros(36, 36) ;

X(2, 4) = -1 ;

X(3, 5) = -2 ;

X(4, 2) = -1 ;
X(4, 23) = +2 ;

X(5, 3) = -1 ;
X(5, 24) = +2 ;

X(7, 19) = +1 ;

X(8, 10) = -1 ;
X(8, 20) = +1 ;

X(9, 21) = +1 ;
X(9, 11) = -2 ;

X(10, 22) = +1 ;
X(10, 29) = +2 ;
X(10, 8) = -1 ;

X(11, 23) = +1 ;
X(11, 30) = +2 ;
X(11, 9) = -1 ;

X(12, 24) = +1 ;

X(13, 25) = +2 ;

X(14, 16) = -1 ;
X(14, 26) = +2 ;

X(15, 17) = -2 ;
X(15, 27) = +2 ;

X(16, 28) = +2 ;
X(16, 35) = +2 ;
X(16, 14) = -1 ;

X(17, 29) = +2 ;
X(17, 36) = +2 ;
X(17, 15) = -1 ;

X(18, 30) = +2 ;

X(19, 7) = +1 ;
X(19, 28) = -2 ;

X(20, 22) = -1 ;
X(20, 29) = -2 ;
X(20, 8) = +1 ;

X(21, 23) = -2 ;
X(21, 30) = -2 ;
X(21, 9) = +1 ;

X(22, 10) = +1 ;
X(22, 20) = -1 ;

X(23, 11) = +1 ;
X(23, 21) = -1 ;

X(24, 12) = +1 ;

X(25, 13) = +1 ;
X(25, 34) = -2 ;

X(26, 28) = -1 ;
X(26, 35) = -2 ;
X(26, 14) = +1 ;

X(27, 29) = -2 ;
X(27, 36) = -2 ;
X(27, 15) = +1 ;

X(28, 26) = -1 ;
X(28, 16) = +1 ;

X(29, 17) = +1 ;
X(29, 27) = -1 ;

X(30, 18) = +1 ;

X(32, 34) = -1 ;

X(33, 35) = -2 ;

X(34, 32) = -1 ;

X(35, 33) = -1 ;

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
   resultC(counter, 1) = result(8, 1) ;        % for ad * a
   
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
