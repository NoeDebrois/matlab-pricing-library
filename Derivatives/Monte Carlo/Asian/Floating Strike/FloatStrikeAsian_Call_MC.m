
% Option considered: Floating Strike Asian Call Option
% Payoff: { S(T) - 1/N * sum_i=1->N ( S(i) )  }⁺

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
DiscPayoff = exp( - r * T ) * max(S(:, end) - mean(S, 2), 0);

% - Compute the price
[price_MC, ~, CI_MC] = normfit( DiscPayoff );


% 2. Antithetic variables technique Monte Carlo ─────────────────────────────────────

% - Simulate prices
[S, Sav] = simulate_BS_AV(S0, T, r, Nsim, N, params);

% - Compute the discounted payoff
DiscPayoff   = exp( - r * T ) * max(S(:, end) - mean(S, 2), 0);
DiscPayoffav = exp( - r * T ) * max(Sav(:, end) - mean(Sav, 2), 0);

% - Compute the price
[price_MCAV, ~, CI_MCAV] = normfit( 0.5 * ( DiscPayoff + DiscPayoffav ) );


% 3. Control variable Monte Carlo ───────────────────────────────────────────────────

% > 3.1. Choose a control variable: f = S(T) - mean(S, 2)

% - expected value of f
Ef = S0 * exp( r * T ) - S0/(N+1) * sum( exp( r * (0:N) * dt ) );

% > 3.2. Estimate alpha

% - Simulate a smaller sample of prices
S = simulate_BS(S0, T, r, Nsim/100, N, params);

f = S(:,end) - mean(S, 2);
g = exp( - r * T ) * max(f, 0);

VC = cov(f, g);

% - Compute alpha
alpha = - VC(1, 2)/VC(1, 1);

% > 3.3. Compute the price

% - Simulate all prices
S = simulate_BS(S0, T, r, Nsim, N, params);

f = S(:,end) - mean(S, 2);
g = exp( - r * T ) * max(f, 0);


% - Compute the price
[price_MCCV, ~, CI_MCCV] = normfit( g + alpha * ( f - Ef ) );


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



