
% Option considered: American Put Option
% Technique: Monte Carlo - Longstaff and Schwartz algorithm
% Payoff: sup_t ∈ [0, T] { K - S(t) }⁺

clc
clear all
close all

%profile on

set(0, 'DefaultFigureWindowStyle', 'docked');


%% Parameters

% option parameters
S0 = 200;      % spot value
r  = 0.1/100;     % risk-free interest rate
T  = 1;        % maturity
K  = 1.2 * S0; % strike

% simulation parameters
Nsim = 1e3;             % number of simulations
%N    = round( T * 365 ); % daily monitoring
%N    = round( T * 252 ); % business day monitoring
N    = round( T * 52 ); % weekly monitoring
%N    = round( T * 12 ); % monthly monitoring
dt   = T/N;             % delta time


%% Calibration

% add path
%addpath('..Your path.../Pricing Library/Models/Black & Scholes')

% model parameters
[params, error_prices, error_vol] = calibrate_BS(S0, r);

SPaths = simulate_BS(S0, T, r, Nsim, N, params);
SPaths = SPaths(:, 2:end); % without S0


%% Initialization

ExerciseTime = N * ones(Nsim, 1);           % initialize the exercise time of the American options with the maturity (= M)
                                         % Note: it could be smaller than the matuirty afterwards
CashFlows    = max(0, K - SPaths(:, N)); % payoff at M (the maturity)


%% Bacward in time procedure

% Note: we consider only the InMoney case, thus option will be exercised

for step = N-1:-1:1 % go backward in time
    
    % only choose in the money options to avoid unnecessary computations
    InMoney = find( K > SPaths(:, step) ); % options in the money at the current step
    S       = SPaths(InMoney, step);       % S(step) of the in the money options
    
    
    % > Regression                    <
    
    % basis functions = [1, S, S^2] having considered M = 3
    RegrMat = [ones(length(S), 1), S, S.^2];                                            % A
    YData   = CashFlows(InMoney) .* exp( - r * dt * ( ExerciseTime(InMoney) - step ) ); % b
    alpha   = RegrMat \ YData;                                                          % x = A \ b
    
    
    % > IV and CV at the current step <
    
    IV = K - S;           % intrinsic value                        <- IV = {K - S}^+
    CV = RegrMat * alpha; % continuation value (its approximation) <- CV = A * x
    
    % > Early Exercise                <
    
    % Paths with early exercise at time step
    Index         = find( IV > CV ); % an option is exercised early if IV > CV since V(step) = max{IV(step), CV(step)}
    ExercisePaths = InMoney(Index);  % index of simulation where it is optimal to early exercise at time step
    
    % Update Cashflows
    CashFlows(ExercisePaths) = IV(Index); % the cashflows of the options exercised early are updated with the IV at the step of exercise
    
    % Update Exercise Time
    ExerciseTime(ExercisePaths) = step; % the exercise times of the options exercised early are updated with the step of exercise
    
end

% compute the price discounting the cashflows
[Price, ~, CI] = normfit( CashFlows .* exp( - r * dt * ExerciseTime ) )
% Remark: the discounted cashflows are the final payoffs if early exercised
% is not used, otherwise are the IV at the step of exercise.
% In the first case the discount is done with as time the time of
% maturity (T = M*dt) otherwise with time equal to step*dt which is the
% time of maturity exercising the option earlier.



