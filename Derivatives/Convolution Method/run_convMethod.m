
clc
clear all
close all

set(0, 'DefaultFigureWindowStyle', 'docked');

% Price a down & out call option using the convolution method

%% Parameters

param.rf    = 0.02; % risk-free rate
param.q     = 0;    % dividend
param.distr = 1;    % distribution chosen (1 for Normal, 2 for NIG)
param.m     = 0;    % drift
param.s     = 0.2;  % volatility
param.T     = 1;    % maturity

S_0      = 1;
K        = 1;
Ndate    = 12;
Barrier  = 0.8;  % lower barrier
N        = 2^12; % parameter for the grid of the Fourier transform

param.dt = param.T/Ndate;


%% Computations

[S, v] = CONV(S_0, K, Ndate, N, Barrier, param);
price  = interp1(S, v, S_0, 'spline')

