function [P_i] = callPrice_FFT_Bates(S0, T, r, K_i, PARAMS)

% Price of a Plain Vanilla Call exploiting the Carr-Madan algorithm
% Model: Bates

%INPUT:
% K_i    : strikes (could be a vector)
%
%OUTPUT:
% P_i    : prices (could be a vector)
%

%% Parameters

% discretization parameters
Npow = 16;
N    = 2^Npow; % grid point
A    = 2000;   % upper bound


%% Computations

% 0. FFT grid ───────────────────────────────────────────────────────────────────────

% integral grid (for dv) to compute integral numerically as a summation
eta  = A/N;
v    = eta*(0:N-1);
v(1) = 1e-22; % correction term
% Note: it cannot be equal to zero, otherwise we get NaN

% logstrike grid, then compute summation via FFT
lambda = 2 * pi/(N * eta); 
k      = - lambda * N/2 + lambda * (0:N-1);
K      = S0 * exp(k); % strike 


% 1. Find the Fourier transform of Z(k) using the Carr-Madan closed formula ─────────

% Fourier transform of z_k
tr_fou = fourier_transform(r, PARAMS, T, v);


% 2. Compute the prices computing the integral using trapezoidal rule ───────────────

% trapezoidal rule
w = [0.5  ones(1, N-2)  0.5];
h = exp( 1i * (0:N-1) * pi ) .* tr_fou .* w * eta;
P = S0 * real( fft(h) / pi + max( 1 - exp( k - r * T ), 0 ) ); % prices


% 3. Compute the call price ─────────────────────────────────────────────────────────

% focus on the grid of prices by deleting too small and too big strikes
index = find( ( K > 0.1 * S0 & K < 3 * S0 ) );
K = K(index);
P = P(index);

% find the option price interpolating on the pricing grid
P_i = interp1(K, P, K_i, 'spline');


end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function fii = fourier_transform(r, PARAMS, T, v)

fii = ( exp( 1i * r * v * T) ) .* ( ( characteristic_func(PARAMS, T, v - 1i) - 1 ) ./ ( 1i * v .* ( 1 + 1i * v ) ) );


end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function f = characteristic_func(PARAMS, T, u)

% Bates characteristic function found as Heston characteristic function plus jumps

% model parameters 
theta  = PARAMS(1); % vol-of-vol
xi     = PARAMS(2); % mean reversion speed
rho    = PARAMS(3); % correlation
eta    = PARAMS(4); % mean 
V0     = PARAMS(5);

mu     = PARAMS(6); % jump mean
delta  = PARAMS(7); % jump standard deviation
lambda = PARAMS(8); % jump intensity
    
alfa     = - 0.5 * ( u .* u + u * 1i );
beta     = xi - rho * theta * u * 1i;
epsilon2 = theta * theta;
gamma    = 0.5 * epsilon2;

D = sqrt(beta .* beta - 4.0 * alfa .* gamma);

bD  = beta - D;
eDt = exp(- D * T);

G   = bD ./ (beta + D);
B   = (bD ./ epsilon2) .* ((1.0 - eDt) ./ (1.0 - G .* eDt));
psi = (G .* eDt - 1.0) ./(G - 1.0);
A   = ((xi * eta) / (epsilon2)) * (bD * T - 2.0 * log(psi));

% jump part
J = @(v) lambda * ( exp( - delta^2 * v.^2/2 + 1i * mu * v ) - 1 ); % characteristic exponent without drift
drift_rn = - J(-1i); % drift chosen under the risk neutral measure

% characteristic exponent
J = drift_rn * 1i * u + J(u);

y = A + B * V0;
f = exp( y + J * T );


end

