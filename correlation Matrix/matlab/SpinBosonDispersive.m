
% **********************************
% *                                *
% * July 23, 2023                   *
% * TED College                    *
% * Spin Boson Model               *
% * 1st Version                    *
% *                                *
% * Here I study the dynamics      *
% * of qubit-three level system    *
% * Dispersive Hamiltonian Version *
% * Here, matrix C is based on     *
% * My own approach                *
% *                                *
% **********************************

tic;
clc;
hold on

% **********************************
% Constants:

alpha = 5;
qDimension = 2;
dimension = 60;

T = 2 * pi ;          % Totoal Time, (2 * pi) for a complete rotation 
step = 1/16 ;          % Powers of 2 : 16, 32, 64, ...

% **********************************
% Variables:

% **************
% Operators: Qubit + Field

% Qubit
vec1 = sqrt(1: qDimension - 1) ;
q = diag(vec1, 1) ;
qd = q' ;
I_qubit = eye(qDimension) ;

% Field
vec2 = sqrt(1: dimension - 1) ;
b = diag(vec2, 1) ;
bd = b' ;
I_oscillator = eye(dimension) ;

sigma = kron(q, I_oscillator) ;
field = kron(I_qubit, b) ;
I = kron(I_qubit, I_oscillator) ;

% **************
% States:

% First of all, let's define our pure states 
ket = zeros(dimension,1) ;
ketZero = ket ;
ketZero(1) = 1;        % |0 > : ground state, [1, 0, 0, ...]
ketOne = ket ;
ketOne(2) = 1 ;        % |1 > : first excited state: [0, 1, 0, ...]
ketTwo = ket ;
ketTwo(3) = 1 ;         % |2 > : second excited state, [0, 0, 1, ...]
ketThree = ket ;
ketThree(4) = 1 ;       % |3 > : third excited state, [0, 0, 0, 1, ...]

qKet = zeros(qDimension, 1) ;
qKetZero = qKet ;
qKetZero(1) = 1 ;       % |0 > : ground state, [1, 0, 0, ...]
qKetOne = qKet ;
qKetOne(2) = 1 ;        % |1 > : first excited state: [0, 1, 0, ...]
qKetTwo = qKet ;
qKetTwo(3) = 1 ;        % |2 > : second excited state, [0, 0, 1, ...]

ketAlpha = exp(-0.5 * abs(alpha) ^ 2) * expm(alpha * bd) * expm(-conj(alpha) * b) * ketZero ;     % coherent state alpha

% Rho = kron(qKetZero * qKetZero', ketAlpha * ketAlpha') ;
Rho = kron(qKetZero * qKetZero', ketOne * ketOne'); 
% Rho = kron(qKetOne * qKetOne', ketTwo * ketTwo');

% **********************************
% Correlation Matrix:

A = cell(1, 3) ;
A_dagger = cell(3, 1) ;
C = cell(3, 3) ;

A{1, 1} = I ;
A{1, 2} = field * (sigma' * sigma) ;
A{1, 3} = field * (I - (sigma' * sigma)) ;

A_dagger{1, 1} = I ;
A_dagger{2, 1} = field' * (sigma' * sigma) ;
A_dagger{3, 1} = field' * (I - (sigma' * sigma)) ;

% =============================================================== change
for i = 1:3
    for j = 1:3
        
        C{i, j} = A_dagger{i, 1} * A{1, j};
    end
end

% **********************************
% Initial value for C:

C_init = zeros(9, 1);
C_final = zeros(9, 1);

k = 0;
for i = 1:3
    for j = 1:3
        
        k = k+1;
        C_init(k) = trace(Rho * C{i,j});
    end
end

% **********************************
% Lindbladian:

W = [0, 0, 0;
     0, 2, 0;
     0, 0, 1];

L = 1i * (kron(W,eye(3)) - kron(eye(3), W.'));

% **********************************
% Main Body:

Time = 0: step: T;

resultC = zeros(round((1/step) * T), 1);

counter = 1;

for time = 0: step: T
    
    co = expm(L * time);
    C_final = co * C_init;
    
    % Populations: 
    resultC(counter, 1) = C_final(5, 1) + C_final(9, 1) ;    % sigma^dagger sigma
    
    disp('--------------------------');
    disp('I am here:');
    disp(C_final(1, 1));
    
    C_final = zeros(9, 1);
    counter = counter + 1;
end

% **********************************
% Plot:

plot(Time, real(resultC));
title('Population vs. Time');
xlabel("Time") ;
ylabel("Population") ;

% **********************************
% Code RunTime:

elapsed_Time = toc;
disp('Run Time (minute):');
disp(elapsed_Time / 60);

