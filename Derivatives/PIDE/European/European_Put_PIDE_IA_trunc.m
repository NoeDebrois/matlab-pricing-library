
clc
clear all
close all

set(0, 'DefaultFigureWindowStyle', 'docked');

% Price an EU Put Option with PIDE method
% Method: Theta Method + Asmussen-Rosinski truncation
% PDE: log-price transform (in x)

%% Parameters

% option parameters
S0 = 200;     % spot value
r  = 0.1/100; % risk-free interest rate
T  = 1;       % maturity
K  = 150;     % strike (ATM option)

% discretization parameters
N = 2000;
M = 200;

% method parameters
theta = 0.5;
epsilon = 0.2;


%% Calibration

% add path
%addpath('..Your path...e/Pricing Library/Models/VG')
%addpath('..Your path.../Pricing Library/Models/Extended VG')
%%addpath('..Your path.../Pricing Library/Models/NIG')
%%addpath('..Your path.../Pricing Library/Models/Extended NIG')

% model parameters
[params, error_prices, error_vol] = calibrate_ExtVG(S0, r);
%[params, error_prices, error_vol] = calibrate_ExtNIG(S0, r);

% significant parameters
sigmaGBM = params(4);


%% Grid

% Note: since we are not dealing with BS model we find the domain
% boundaries approximatively, not using an approximation of the normal distribution

% space
x_min = log( 0.2 * S0 / S0 ); % S_min = 0.2 * S0
x_max = log(3);               % S_max = 3 * S0
x = linspace(x_min, x_max, N+1);
dx = x(2) - x(1);

% time
dt = T/M;


%% Compute alpha and lambda

% Levy measure
nu = @(y) LevyMeasure_VG(y, params);
%nu = @(y) LevyMeasure_NIG(y, params);

% integrals that can be computed if working with finite acitivty levy processes
[alpha, lambda, sigma_eps, LB, UB] = levy_integral(nu, N, epsilon);

% truncated Levy measure
nu = @(y) nu(y) .* ( abs(y) > epsilon );

% plot
figure
y = linspace(LB, UB, N);
plot(y, nu(y));


%% Matrix

% Note: we use the same coefficients as in the Black-Scholes with the added lambda and alpha since
% we approximate the process and take it as if it was a finite activity process

% model sigma
sigma = sqrt( sigmaGBM^2 + sigma_eps^2 );

% coefficients
A = (1 - theta) * ( - ( r - sigma^2/2 - alpha )/( 2 * dx ) + sigma^2/( 2 * dx^2 ) );
B = - 1/dt + (1 - theta) * ( - sigma^2/( dx^2 ) - ( r + lambda ) );
C = + (1 - theta) * ( ( r - sigma^2/2 - alpha )/( 2 * dx ) + sigma^2/( 2 * dx^2 ) );

At = - theta * ( - ( r - sigma^2/2 - alpha )/( 2 * dx ) + sigma^2/( 2 * dx^2 ) );
Bt = - 1/dt - theta * ( - sigma^2/( dx^2 ) - ( r + lambda ) );
Ct = - theta * ( ( r - sigma^2/2 - alpha )/( 2 * dx ) + sigma^2/( 2 * dx^2 ) );

% build the matrix
M1 = sparse(N+1, N+1);
M2 = sparse(N+1, N+1);
M1(1, 1) = 1;
M1(end, end) = 1;
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

    % Compute the Integral at time t_j
    I = levy_integral2(LB, UB, x, V, nu, S0, K, exp( - r * ( T - j * dt ) ));

    % Solve the linear system
    V = M1 \ ( M2 * V + BC - I );

end


%% Output

% plot
figure
set(gcf, 'Color', 'w')

plot(x, V);

title('Solution');

% price
price = interp1(x, V, 0, 'spline'); % log(S0/S0) = 0

% Carr & Madan price
CM_price = callPrice_FFT_ExtVG(S0, T, r, K, params);
%CM_price = callPrice_FFT_ExtNIG(S0, T, r, K, params);

% put-call pairity to find put price
CM_price = CM_price - S0 + K * exp( - r * T );

% error
error = CM_price - price;

fprintf('\nPrice               : %.5f\n', price)
fprintf('Carr & Madan price  : %.5f\n', CM_price)
fprintf('Error               : %.5f\n\n', error)


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function [alpha, lambda_num, sigma_eps, LB, UB] = levy_integral(nu, N, epsilon)

% Compute the "known" integrals that sum up to the linear part

% 1. Integral domain truncation
step = 0.5;
tol = 10^-10;

% find the lower bound as the value such that the levy measure is bigger than a tolerance
LB = - max(step, epsilon);
% Note: take the maximum between step and epsilon since before epsilon nu is equal to zero
while nu(LB) > tol
    LB = LB - step;
end

% find the upper bound as the value such that the levy measure is bigger than a tolerance
UB = + max(step, epsilon);
while nu(UB) > tol
    UB = UB + step;
end

% 2. Compute the integral using quadrature
N_q = N; % number of quadrature points (of the laft and the right intervals separately)
% Note: when we improve the accuracy on the number of steps the domain is divided into
% we improve the accuracy of the numerical quadrature

y = linspace(LB, - epsilon, N_q); % quadrature left domain
alpha = trapz(y, ( exp(y) - 1 ) .* nu(y)); % left alpha
lambda_num = trapz(y, nu(y)); % left lambda

% plot left
figure
plot(y, nu(y))
title('Levy measure')

y = linspace(epsilon, UB, N_q); % quadrature right domain
alpha = alpha + trapz(y, ( exp(y) - 1 ) .* nu(y)); % left+right alpha
lambda_num = lambda_num + trapz(y, nu(y)); % left+right lambda

% plot right
hold on
plot(y, nu(y))

% sigma(epsilon)
y = linspace( - epsilon, epsilon, 2 * N ); % N even
sigma_eps = sqrt( trapz(y, y.^2 .* nu(y) ) );

end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function I = levy_integral2(LB, UB, x, V, nu, S0, K, disc)

% Compute the "unknown" integral that constitutes the non linear part

% initialization
I = zeros(size(V));
% Note: in this way I(1) = I(end) = 0 to guarantee to have the boundary conditions

% domain
N_q  = length(x) - 1; % number of quadrature points
% Note: we take an even number of quadrature points since we are working
% with a symmetric domain to avoid having nu(0) which is equal to Inf with
% infinite activity processes 
y    = linspace(LB, UB, N_q)';
dy   = y(2) - y(1);
nu_y = nu(y);

% trapezoidal quadrature
w      = ones(N_q, 1); % weights
w(1)   = 0.5;
w(end) = 0.5;

for i = 2:length(I)-1 % it stays that I(1) = I(end) = 0

    I(i) = sum( w .* V_f(x, V, x(i) + y, S0, K, disc) .* nu_y ) * dy;

end


end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function v = V_f(x, V, y, S0, K, disc)

% Function that gives the value function in the point of interest

% initialization
v = zeros(size(y));

% 1. y <= x_min
index    = ( y <= x(1) );
v(index) = K * disc - S0 * exp( y(index) ); % BC for y <= x_min

% 2. y >= x_max
index    = ( y >= x(end) );
v(index) = 0; % BC for y >= x_max

% 3. x_min < y < x_max
index    = find( ( y > x(1) ) .* ( y < x(end) ) );
v(index) = interp1(x, V, y(index)); % value for x_min < y < x_max


end

