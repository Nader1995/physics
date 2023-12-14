
% **********************************
% *                                *
% * June 25, 2023                  *
% * TED College                    *
% * Spin Boson Model               *
% * 1st Version                    *
% *                                *
% * Here I study the dynamics      *
% * of qubit-two level system      *
% *                                *
% **********************************

tic;
clc;
hold on

% **********************************
% Constants:

d_qubit = 2;
d_oscillator = 4;

T = 2 * pi ;          % Totoal Time, (2 * pi) for a complete rotation 
step = 1/8 ;          % Powers of 2 : 16, 32, 64, ...

% **********************************
% Variables:

% **************
% States:

psi_qubit = zeros(d_qubit, 1);
psi_qubit(1) = 1;

psi_oscillator = zeros(d_oscillator, 1);
psi_oscillator(4) = 1;

psi = kron(psi_qubit, psi_oscillator);
Rho = psi * psi';
% psi'

% **************
% Operators:

vec1 = sqrt(1: d_qubit - 1);
sigma = diag(vec1, +1);
I_qubit = eye(d_qubit);
% sigma'

vec2 = sqrt(1: d_oscillator - 1);
a = diag(vec2, +1);
I_oscillator = eye(d_oscillator);
% a'

% **********************************
% Correlation Matrix:

A = cell(1, 4);
C = cell(4, 4);

A{1, 1} = kron(I_qubit, I_oscillator);
A{1, 2} = kron(I_qubit, a);
A{1, 3} = kron(sigma, I_oscillator);
A{1, 4} = kron(sigma, a);

for i = 1:4
    for j = 1:4
        
        C{i, j} = A{1, i}' * A{1, j};
    end
end

% **********************************
% Initial value for C:

C_init = zeros(16, 1);

k = 0;
for i = 1:4
    for j = 1:4
        
        k = k+1;
        C_init(k) = trace(Rho * C{i,j});
    end
end

% **********************************
% Lindbladian:
 
W = [0, 0, 0, 0;
     0, 0, 1, 0;
     0, 1, 0, 0;
     0, 0, 0, 0];

L = 1i * (kron(W,eye(4)) - kron(eye(4), W.'));
% **********************************
% Main Body:

Time = 0: step: T;

% Stores the population of photon inside the oscillator
resultC_a = zeros(round((1/step) * T), 1);

% Stores the population of photon inside the qubit
resultC_sigma = zeros(round((1/step) * T), 1);

counter = 1;

for time = 0: step: T
    
    co = expm(L * time);
    result = co * C_init;
    
    % Populations: 
    resultC_a(counter, 1) = result(6, 1);     % a^dagger a
    resultC_sigma(counter, 1) = result(11, 1);    % sigma^dagger sigma
    
    disp('--------------------------');
    disp('I am here:');
    disp(time);
    
    counter = counter + 1;
end
% disp(resultC_a);
% disp(resultC_sigma);

% **********************************
% Plot:

plot(Time, real(resultC_a));
title('Population vs. Time');
xlabel("Time") ;
ylabel("Population") ;

plot(Time, real(resultC_sigma));
title('Population vs. Time');
xlabel("Time") ;
ylabel("Population") ;

% **********************************
% Code RunTime:

elapsed_Time = toc;
disp('Run Time (minute):');
disp(elapsed_Time / 60);
