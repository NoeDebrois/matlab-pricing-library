
clc
clear all
close all

set(0, 'DefaultFigureWindowStyle', 'docked');

% Price a Up and Out Call Option in the Black & Scholes framework
% Method: Theta Method
% PDE: log-price transform (in x)

%% Parameters

% option parameters
S0 = 199.80;   % spot value
r  = 0.1/100;  % risk-free interest rate
T  = 1;        % maturity
K  = S0;       % strike (ATM option)
U  = 1.3 * S0; % barrier

% discretization parameters
N = 2000;
M = 100;

% method parameters
theta = 0.5;


%% Calibration

% add path
%addpath('..Your path.../Pricing Library/Models/Black & Scholes')

% model parameters
[params, error_prices, error_vol] = calibrate_BS(S0, r);

sigma = params;


%% Grid

% log barrier
u = log( U / S0 );

% space
x_min = log( 0.2 * S0 / S0 ); % S_min = 0.2 * S0
x_max = u;
x = linspace(x_min, x_max, N+1);
dx = x(2) - x(1);

% time
dt = T/M;


%% Matrix

% coefficients
A = (1 - theta) * ( - ( r - sigma^2/2 )/( 2 * dx ) + sigma^2/( 2 * dx^2 ) );
B = - 1/dt + (1 - theta) * ( - sigma^2/( dx^2 ) - r );
C = + (1 - theta) * ( ( r - sigma^2/2 )/( 2 * dx ) + sigma^2/( 2 * dx^2 ) );

At = - theta * ( - ( r - sigma^2/2 )/( 2 * dx ) + sigma^2/( 2 * dx^2 ) );
Bt = - 1/dt - theta * ( - sigma^2/( dx^2 ) - r );
Ct = - theta * ( ( r - sigma^2/2 )/( 2 * dx ) + sigma^2/( 2 * dx^2 ) );

% build the matrix
M1 = sparse(N+1, N+1);
M2 = sparse(N+1, N+1);
M1(1, 1) = 1;
M1(end, end)=1;
for i = 1:N-1
     M1(i+1, [i i+1 i+2]) = [A B C];
     M2(i+1, [i i+1 i+2]) = [At Bt Ct];
end


%% Backward in time solution

% at maturity
V = max(S0 * exp(x') - K, 0) .* ( x' < u );

% initialize boundary condition vector
BC = zeros(N+1, 1);

for j = M:-1:1 % known t_j --> unknown t_{j-1}

    % Remark: BC is always zero for an Up & Out call option

    % solve the system
    V = M1 \ ( M2 * V + BC );

end


%% Output

% plot
figure
set(gcf, 'Color', 'w')

plot(x, V);

title('Solution');

% price
price = interp1(x, V, 0, 'spline'); % log(S0/S0) = 0

fprintf('\nPrice        : %.5f\n', price)



