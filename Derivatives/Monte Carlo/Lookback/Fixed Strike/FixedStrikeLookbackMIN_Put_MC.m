
% Option considered: Fixed Strike Lookback Put on Minimum
% Payoff: { K - m(T) }⁺ where m(T) = min_t ∈ [0, T] ( S(t) )

clc
clear all
close all

%profile on

set(0, 'DefaultFigureWindowStyle', 'docked');


%% Parameters

% option parameters
S0 = 199.80; % spot value
r  = 0.01;   % risk-free interest rate
T  = 1;      % maturity
K  = S0;     % strike 

% simulation parameters
Nsim = 1e6;              % number of simulations
%N    = round( T * 365 ); % daily monitoring
%N    = round( T * 252 ); % business day monitoring
N    = round( T * 52 );  % weekly monitoring
%N    = round( T * 12 );  % monthly monitoring
dt   = T/N;              % delta time


%% Calibration

% add path
%addpath('..Your path.../Pricing Library/Models/Black & Scholes')

% model parameters
[params, error_prices, error_vol] = calibrate_BS(S0, r);


%% Pricing

% 1. Classical Monte Carlo ──────────────────────────────────────────────────────────

% - Simulate prices
S = simulate_BS(S0, T, r, Nsim, N, params);

% - Compute the discounted payoff
DiscPayoff = exp( - r * T ) * max(K - min(S, [], 2), 0);

% - Compute the price
[price_MC, ~, CI_MC] = normfit( DiscPayoff );


% 2. Antithetic variables technique Monte Carlo ─────────────────────────────────────

% - Simulate prices
[S, Sav] = simulate_BS_AV(S0, T, r, Nsim, N, params);

% - Compute the discounted payoff
DiscPayoff   = exp( - r * T ) * max(K - min(S, [], 2), 0);
DiscPayoffav = exp( - r * T ) * max(K - min(Sav, [], 2), 0);

% - Compute the price
[price_MCAV, ~, CI_MCAV] = normfit( 0.5 * ( DiscPayoff + DiscPayoffav ) );


% 3. Control variable Monte Carlo ───────────────────────────────────────────────────

% GOOD CONTROL VARIABLE NOT FOUND

% - Compute the price
price_MCCV = zeros(size(price_MCAV));
CI_MCCV    = zeros(size(CI_MCAV));


%% Check Risk-Neutrality

[check, ~, check_IC] = normfit( exp( - r * T ) * S(:, end) );
[check_AV, ~, check_IC_AV] = normfit( exp( - r * T ) * Sav(:, end) );


fprintf('\nCheck risk-neutrality :\n\n')
fprintf('              S0     -        Approximation \n')
fprintf('  MC     :  %.2f   -   %.3f  in  (%.3f, %.3f)\n', S0, check, check_IC(1), check_IC(2))
fprintf('  MC AV  :  %.2f   -   %.3f  in  (%.3f, %.3f)\n\n', S0, check_AV, check_IC_AV(1), check_IC_AV(2))


%% Output

fprintf('\n')
disp('         ┌───── Price ──────────── CI ─────────── |CI| ─────┐')
fprintf('  MC     │     %.5f     (%.5f, %.5f)    %.5f    │\n', price_MC, CI_MC(1), CI_MC(2), CI_MC(2)-CI_MC(1))
fprintf('  MC AV  │     %.5f     (%.5f, %.5f)    %.5f    │\n', price_MCAV, CI_MCAV(1), CI_MCAV(2), CI_MCAV(2)-CI_MCAV(1))
fprintf('  MC CV  │     %.5f     (%.5f, %.5f)    %.5f    │\n', price_MCCV, CI_MCCV(1), CI_MCCV(2), CI_MCCV(2)-CI_MCCV(1))
disp('         └──────────────────────────────────────────────────┘')
fprintf('\n')

% profiler
%profile viewer
%profile off



