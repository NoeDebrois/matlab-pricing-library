function PRICE = callPrice_FFT_Kou(S0, T, r, STRIKE, PARAMS)

% Price of a Plain Vanilla Call exploiting the Carr-Madan algorithm
% Model: Kou

format long


%% Parameters

% model parameters
SIGMA        = PARAMS(1); % standard deviation of the Brownian motion
LAMBDA       = PARAMS(2); % intensity of the jumps
LAMBDA_MINUS = PARAMS(3); % 1/mean of the positive jump size
LAMBDA_PLUS  = PARAMS(4); % 1/mean of the negative jump size
P            = PARAMS(5); % probability of having a positive jump

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
CharFunc = @(v) exp( T * CharExp_Kou(v, SIGMA, LAMBDA, LAMBDA_MINUS, LAMBDA_PLUS, P) );

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


function V = CharExp_Kou(v, SIGMA, LAMBDA, LAMBDA_MINUS, LAMBDA_PLUS, P)

% risk-neutral characteristic exponent

V = @(v) - SIGMA^2 .* v.^2/2 + 1i .* v .* LAMBDA .* ( P ./ ( LAMBDA_PLUS - 1i .* v ) - ( 1 - P ) ./ ( LAMBDA_MINUS + 1i .* v ) ); % characteristic exponent without drift
drift_rn = - V(-1i); % drift chosen under the risk neutral measure (the char exp is equal to zero)

% characteristic exponent
V = drift_rn * 1i * v + V(v);


end

