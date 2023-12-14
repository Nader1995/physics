
% **********************************
% *                                *
% * June 8, 2023                   *
% * TED College                    *
% * Spin Boson Model               *
% * 2st Version                    *
% *                                *
% * Here I study the dynamics      *
% * of qubit-three level system    *
% *                                *
% *                                *
% **********************************

tic;
clc;
hold on

% **********************************
% Constants:

d_qubit = 2;

% =============================================================== change
d_oscillator = 3;

T = 2 * pi ;          % Totoal Time, (2 * pi) for a complete rotation 
step = 1/8 ;          % Powers of 2 : 16, 32, 64, ...

% **********************************
% Variables:

% **************
% States:

psi_qubit = zeros(d_qubit, 1);
psi_qubit(1) = 1;

psi_oscillator = zeros(d_oscillator, 1);
% =============================================================== change
psi_oscillator(3) = 1;

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

% =============================================================== change
A = cell(1, 6);
C = cell(6, 6);

% =============================================================== change
A{1, 1} = kron(I_qubit, I_oscillator);
A{1, 2} = kron(I_qubit, a);
A{1, 3} = kron(I_qubit, a * a);
A{1, 4} = kron(sigma, I_oscillator);
A{1, 5} = kron(sigma, a);
A{1, 6} = kron(sigma, a * a); 

% =============================================================== change
for i = 1:6
    for j = 1:6
        
        C{i, j} = A{1, i}' * A{1, j};
    end
end

% **********************************
% Initial value for C:

% =============================================================== change
C_init = zeros(36, 1);

% =============================================================== change
k = 0;
for i = 1:6
    for j = 1:6
        
        k = k+1;
        C_init(k) = trace(Rho * C{i,j});
    end
end

% **********************************
% Lindbladian:

% =============================================================== change
W = [0, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 0, 0];

% =============================================================== change
L = 1i * (kron(W,eye(6)) - kron(eye(6), W.'));

% **********************************
% Main Body:

Time = 0: step: T;

resultC_a = zeros(round((1/step) * T), 1);
resultC_sigma = zeros(round((1/step) * T), 1);

counter = 1;

for time = 0: step: T
    
    co = expm(L * time);
    result = co * C_init;
    
    % =============================================================== change
    % Populations: 
    resultC_a(counter, 1) = result(8, 1);     % a^dagger a
    resultC_sigma(counter, 1) = result(22, 1);    % sigma^dagger sigma
    
    disp('--------------------------');
    disp('I am here:');
    disp(time);
    
    counter = counter + 1;
end
% disp(resultC_a);
% disp(resultC_sigma);

% **********************************
% Plot:

% plot(Time, real(resultC_a));
% title('Population vs. Time');
% xlabel("Time") ;
% ylabel("Population") ;

plot(Time, real(resultC_sigma));
title('Population vs. Time');
xlabel("Time") ;
ylabel("Population") ;

% **********************************
% Code RunTime:

elapsed_Time = toc;
disp('Run Time (minute):');
disp(elapsed_Time / 60);
