function [X, Xav] = simulate_BS_AV(S0, T, r, N_SIM, N, PARAMS)

% Simulate Black-Scholes prices with antithetic variables Monte Carlo technique starting from the logreturns X
% S(t) = S0 * exp( rt + X(t) )

%% Parameters

% model parameters
SIGMA = PARAMS;

% delta time
dt = T/N;


%% Computations

% initialize
X   = zeros(N_SIM, N+1); % standard Monte Carlo
Xav = X;                 % antithetic variables Monte Carlo

% sample standard normal
Z = randn(N_SIM, N); % sample all random variables
% Note: if Nsim and N are too big we could have an out of memory issue
% since it would impossible to store such a big matrix

for i = 1:N

    X(:, i+1)   = X(:, i) + ( r - SIGMA^2/2 ) * dt + SIGMA * sqrt(dt) * Z(:, i);
    Xav(:, i+1) = Xav(:, i) + ( r - SIGMA^2/2 ) * dt - SIGMA * sqrt(dt) * Z(:, i);

end

% from logreturns to prices: X -> S
X   = S0 * exp( X );
Xav = S0 * exp( Xav );


end

