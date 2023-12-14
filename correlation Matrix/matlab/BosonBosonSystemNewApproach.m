
% **********************************
% *                                *
% * July 12, 2023                  *
% * TED College                    *
% * Boson-Boson Model,             *
% * The new approach 1st Version   *
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
dimension2 = 2 ;       % Dimension of restricted bosonic subspace  

% =======================================================================
% Defining Annihilation and Creation Operators 

% Bosonic Subspace
bI = eye(dimension1) ;
bVec = sqrt(1:dimension1 - 1) ;
b = diag(bVec, 1) ;
bd = b' ;

% Atomic Subspace
aI = eye(dimension2) ;
aVec = sqrt(1:dimension2 - 1) ;
a = diag(aVec, 1) ;
ad = a' ;

% =======================================================================
% Defining Annihilation and Creation Operators of Bosonic and Atomic Fields
% The Order of Subspaces is as follows:
% Atomic + Bosonic = A + B

A = kron(a, bI) ;
B = kron(aI, b) ;

AD = A' ;
BD = B' ;

% =======================================================================
% Defining the States

% First of all, let's define our pure states 
bket = zeros(dimension1,1) ;
bketZero = bket ;
bketZero(1) = 1;        % |0 > : ground state, [1, 0, 0, ...]
bketOne = bket ;
bketOne(2) = 1 ;        % |1 > : first excited state: [0, 1, 0, ...]
bketTwo = bket ;
bketTwo(3) = 1 ;         % |2 > : ground state, [0, 0, 1, ...]
bketThree = bket ;
bketThree(4) = 1 ;       % |3 > : ground state, [0, 0, 0, 1, ...]

aket = zeros(dimension2, 1) ;
aketZero = aket ;
aketZero(1) = 1 ;       % |0 > : ground state, [1, 0, 0, ...]
aketOne = aket ;
aketOne(2) = 1 ;
 
ketAlpha = exp(-0.5 * abs(alpha) ^ 2) * expm(alpha * bd) * expm(-conj(alpha) * b) * bketZero ;     % coherent state alpha

% state = kron(bketZero * bketZero', ketAlpha * ketAlpha') ;
psi = kron(aketZero, bketThree) ;
state = kron(aketZero * aketZero', bketTwo * bketTwo') ; 

C = cell(2, 2) ;
res = cell(2, 2) ;
finalResult = cell(1, 1) ;
finalResult{1, 1} = zeros(120) ;

C{1, 1} = AD * A ;
C{1, 2} = AD * B ;
C{2, 1} = BD * A ;
C{2, 2} = BD * B ;

C_init = zeros(4, 1) ;

k = 0 ;
for i = 1:2
    for j = 1:2
      
    res{i, j} = zeros(120) ;
    k = k + 1 ;
    C_init(k) = trace(state * C{i, j}) ;
    
    end
end

W = [ 0, 1 ;
      1, 0 ] ;

lindblad = -1i * (kron(W, eye(2)) - kron(eye(2), W.')) ;

% Dynamics

Time = 0: step: T ;
resultC = zeros(round((1/step) * T), 1) ;

counter = 1 ;

% for time = 0: step: T
%    
%    co = expm(lindblad * time);
%    result = co * C_init;
%     
%    resultC(counter, 1) = result(1, 1);
%    
%    disp( " ================================================ " );
%    disp( resultC(counter, 1) );
%  
%    disp("I am here: ");
%    disp(time);
%     
%    counter = counter + 1;
%     
% end

for time = 0: step: T
    
    % ====================================
    % Time evolution: Distinguishable case
    U = expm(-1i * W * time);
    UDag = expm(1i * W' * time);
    
    for i = 1 : 2
        for k = 1 : 2
            
            % first: U * C
            res{1, i} = res{1, i} + U(1, k) * C{k, i};            
        end
    end
    
    for i = 1 : 2
       
        % second: C * UDag
        finalResult{1, 1} = finalResult{1, 1} + res{1, i} * UDag(i, 1);
    end
   
    resultC(counter, 1) = psi' * finalResult{1, 1} * psi;
    
    counter = counter + 1;
    
    for i = 1 : 2
        for k = 1 : 2
        
        res{i, k} = zeros(120);
        end
    end
    
    finalResult{1, 1} = zeros(120);
    
end

% =======================================================================
% Plot the Result                         

plot(Time, real(resultC));
title('Population ')
xlabel("Time") ;
ylabel("<a dagger a>") ;

elapsed_Time = toc ;
disp('Run Iime (minute):') ;
disp(elapsed_Time / 60);




