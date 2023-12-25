
clc
clear all
close all

set(0, 'DefaultFigureWindowStyle', 'docked');

% Price an EU Put Option in the Black & Scholes framework
% Method: Theta Method
% PDE: log-price transform (in x)

%% Parameters

% option parameters
S0 = 199.80;   % spot value
r  = 0.1/100;  % risk-free interest rate
T  = 1;        % maturity
K  = S0;       % strike (ATM option)

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

% space
x_min = ( r - sigma^2/2 ) * T - 6 * sigma * sqrt(T);
x_max = ( r - sigma^2/2 ) * T + 6 * sigma * sqrt(T);
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

% value of the option at maturity
V = max(K - S0 * exp(x'), 0);

% initialize boundary condition vector
BC = zeros(N+1, 1);

for j = M:-1:1 % known t_j --> unknown t_{j-1}

    % boundary condition
    BC(end) = K * exp( - r * ( T - (j - 1) * dt ) ) - S0 * exp(x_min);

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

% exact price
ex_price = blsprice(S0, K, r, T, sigma);

% error
error = ex_price - price;

fprintf('\nPrice        : %.5f\n', price)
fprintf('Exact price  : %.5f\n', ex_price)
fprintf('Error        : %.5f\n\n', error)



