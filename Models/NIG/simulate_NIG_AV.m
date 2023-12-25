function [S, Sav] = simulate_NIG_AV(S0, T, r, N_SIM, N, PARAMS)

% Simulate Normal Inverse Gaussian prices with Monte Carlo antithetic variables technique
% starting from Inverse Gaussian subordinators
% S(t) = S0 * exp( rt + X(t) )
% X(t) is a Normal Inverse Gaussian process in the Risk Neutral measure

%% Parameters

% model parameters
SIGMA = PARAMS(1);
THETA = PARAMS(2);
K_NIG = PARAMS(3);

% delta time
dt = T/N;


%% Computations

% > inverse Gaussian simulation

% 1. icdf simulation
%G = icdf('InverseGaussian', rand(N_SIM, N), dt, dt^2/K_NIG);

% 2. Analogic simulation

% generate independent uniform and chi-squared random variables
U = unifrnd(0, 1, [N_SIM, N]); % UNIF(0, 1)
Y = randn(N_SIM, N) .^ 2;      % chi-squared

% find Gstar from u and x2
mu = dt;
lambda = dt^2 / K_NIG;
X1 = mu + 0.5 .* mu.^2 .* Y ./ lambda - 0.5 * mu ./ lambda * sqrt( 4 .* mu .* lambda .* Y + mu.^2 * Y .^2 );
temp = mu ./ ( X1 + mu );
indexes1 = U <= temp;
indexes2 = U > temp;
% Note: indexes1 and indexes2 are matrices of indexes

% build G
G           = zeros(N_SIM, N);
G(indexes1) = X1(indexes1);
G(indexes2) = mu.^2 ./ X1(indexes2);


% > characteristic exponent

charexp = @(u) 1 / K_NIG - 1 / K_NIG * sqrt( 1 + u^2 * SIGMA^2 * K_NIG - 2i * THETA * K_NIG * u );
drift = r - charexp(-1i); % in order to be risk-neutral


% > compute logreturns

% initialization
X = zeros(N_SIM, N+1);
Xav = zeros(N_SIM, N+1);

% sample standard normal
Z = randn(N_SIM, N);

for i = 1:N
    
    X(:, i+1)   = X(:, i)   + drift * dt + THETA * G(:, i) + SIGMA * sqrt(G(:, i)) .* Z(:, i);
    Xav(:, i+1) = Xav(:, i) + drift * dt + THETA * G(:, i) - SIGMA * sqrt(G(:, i)) .* Z(:, i);

end

S = S0 * exp( X );
Sav = S0 * exp( Xav );


end

