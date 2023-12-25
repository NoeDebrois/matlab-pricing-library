
clc
clear all
close all

set(0, 'DefaultFigureWindowStyle', 'docked');

% Price an EU Call Option in the Black & Scholes framework
% Method: Implicit Euler
% PDE: classic (in S)
% Solve the linear system with SOR

%% Parameters

% option parameters
S0 = 200;     % spot value
r  = 0.1/100; % risk-free interest rate
T  = 1;       % maturity
K  = 150;     % strike (ATM option)

% discretization parameters
N = 500;
M = 200;

% method parameters
theta = 0.5;

% SOR parameters
maxiter = 500;
tol     = 1e-4;
w       = 1.5;


%% Calibration

% add path
%addpath('..Your path.../Pricing Library/Models/Black & Scholes')

% model parameters
[params, error_prices, error_vol] = calibrate_BS(S0, r);

sigma = params;


%% Grid

% space
% S_min = S0 * exp( ( r - sigma^2/2 ) * T - 6 * sigma * sqrt(T) );
% S_max = S0 * exp( ( r - sigma^2/2 ) * T + 6 * sigma * sqrt(T) );
S_min = 0.1 * S0; % smaller truncation
S_max = 5 * S0;   % smaller truncation
S = linspace(S_min, S_max, N+1);
dS = S(2) - S(1);

% time
dt = T/M;


%% Matrix

% S_i values
node = S(2:end-1)';

% coefficients
A = - r * node/( 2 * dS ) + ( sigma * node ).^2/( 2 * dS^2 );
B = - 1/dt - ( sigma * node ).^2/( dS^2 ) - r;
C = + r * node/( 2 * dS ) + ( sigma * node ).^2/( 2 * dS^2 );

% system matrix
Mat = sparse(N+1, N+1);
Mat(1, 1) = 1;
for i = 1:N-1
    Mat(i+1, [i i+1 i+2]) = [A(i) B(i) C(i)];
end
Mat(N+1, N+1) = 1;


%% Backward in time solution

% at maturity
V = max(K - S', 0);

for j = M:-1:1 % known t_j --> unknown t_{j-1}

    rhs = [ K - S_min;
            - 1/dt * V(2:end-1);
            0 ];

    Vold = V; % initial guess
    for k = 1:maxiter

        for i = 1:N+1

            if i == 1
                y = 1/Mat(i, i) * ( rhs(i) - Mat(i, i+1) * Vold(i+1) );
            elseif i == N+1
                y = 1/Mat(i, i) * ( rhs(i) - Mat(i, i-1) * V(i-1) );
            else
                y = 1/Mat(i, i) * ( rhs(i) - Mat(i, i+1) * Vold(i+1) - Mat(i, i-1) * V(i-1) );
            end
            V(i) = max( Vold(i) + w * ( y - Vold(i) ), K - S(i) );
            % V(t, S) >= K - S inside the linear system
        end
        if norm(V - Vold, 'inf') < tol
            [j k]
            break;
        else
            Vold = V;
        end
   end

end


%% Output

% plot
figure
set(gcf, 'Color', 'w')

plot(S,V);

title('Solution');

% price
price = interp1(S, V, S0, 'spline'); % log(S0/S0) = 0

fprintf('\n\nPrice        : %.5f\n', price)


%% Computing greeks

Delta = ( V(3:end) - V(1:end-2) ) / ( 2 * dS );
Gamma = ( V(3:end) - 2 * V(2:end-1) + V(1:end-2) ) / ( dS^2 );

% plot
figure
set(gcf, 'Color', 'w')

subplot(1, 2, 1)
hold on

plot(node, Delta)
% Note: use node since we need N-1 values of S (which are i = 1, ..., N-1, MATLAB i = 2, ..., N)

title('Delta')

subplot(1, 2, 2)
hold on

plot(node, Gamma)

title('Gamma')



