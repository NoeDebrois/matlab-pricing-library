function PRICE = callPrice_FFT_ExtNIG(S0, T, r, STRIKE, PARAMS)

% Price of a plain vanilla call exploiting the Carr-Madan algorithm
% Model: Extended Normal Inverse Gaussian

format long


%% Parameters

% model parameters
SIGMA 	  = PARAMS(1); % standard deviation of the subordinated Brownian Motion
THETA     = PARAMS(2); % drift of the subordinated Brownian Motion
K_NIG     = PARAMS(3); % variance of the subordinator
SIGMA_BAR = PARAMS(4); % standard deviation of the added Brownian Motion

% discretization parameters
Npow = 16;
N    = 2^Npow; % discretization points
A    = 2000;   % domain upper boundary


%% Computations

% 0. FFT grid ───────────────────────────────────────────────────────────────────────

% integral domain grid (for dv)... to compute integral numerically as a
% summation with a quadrature rule (trapezoidal)
eta  = A/N;               % spacing on the domain grid
v    = [0:eta:A*(N-1)/N]; % integral domain grid (only positive half)
v(1) = 1e-22;             % adjust starting point near zero
% Note: the correction is necessary, v(1) could not be equal to zero
% otherwise we get a NaN

% logstrike grid (for dK)... to then compute summation via FFT
lambda = 2 * pi/(N * eta);                  % coefficient
k      = - lambda * N/2 + lambda * (0:N-1); % logstrikes grid


% 1. Find the Fourier transform of Z(k) using the Carr-Madan closed formula ─────────

% model characteristic function
CharFunc = @(v) exp( T * CharExp_ExtNIG(v, SIGMA, THETA, K_NIG, SIGMA_BAR) );

% Fourier transform of Z(k)
g = exp( 1i * r * v * T ) .* ( CharFunc(v - 1i) - 1 )./( 1i * v .* ( 1 + 1i * v ) );


% 2. Compute Z(k) applying the anti Fourier transform ───────────────────────────────

% compute Z(k) on the grid
w = ones(1, N); w(1) = 0.5; w(end) = 0.5;      % trapezoidal rule weights
x = w .* eta .* g .* exp( 1i * pi * (0:N-1) ); % function in the DFT
z_k = real( fft(x) / pi );                     % apply discrete Fourier transform
% Note: division by pi since we are considering only the positive half of
% the domain since Z(k) must be real


% 3. Compute the call price exploiting the value of Z(k) ────────────────────────────

C = S0 * (z_k + max(1 - exp(k - r * T), 0)); % call prices
K = S0 * exp(k);                             % strikes

% focus on the grid of prices (delete too small and too big strikes)
index = find( K > 0.1 * S0 & K < 3 * S0 );
C = C(index);
K = K(index);

% find the option price interpolating on the pricing grid in the strike
% given on the set of strikes K present on the grid
PRICE = interp1(K, C, STRIKE, 'spline');


end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function V = CharExp_ExtNIG(u, SIGMA, THETA, K_NIG, SIGMA_BAR)

% characteristic exponent without drift
V = @(v) - SIGMA_BAR^2/2 * v.^2 + 1 / K_NIG - 1 / K_NIG * sqrt( 1 + v.^2 * SIGMA.^2 * K_NIG - 2i * THETA * K_NIG * v );

% drift chosen under the risk neutral measure (the char exp is equal to zero)
drift_rn = - V(-1i);

% characteristic exponent
V = drift_rn * 1i * u + V(u);


end

