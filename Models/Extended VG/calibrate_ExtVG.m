function [PARAMS, ERROR_PRICES, ERROR_VOL] = calibrate_ExtVG(S0, r)

% Calibrate Extended VG Model using Carr-Madan algorithm for pricing

% SIGMA     = PARAMS(1) - standard deviation of the subordinated Brownian Motion
% THETA     = PARAMS(2) - drift of the subordinated Brownian Motion
% K_VG      = PARAMS(3) - variance of the subordinator
% SIGMA_BAR = PARAMS(4) - standard deviation of the added Brownian Motion


%% Input

% add paths
%addpath('..Your path.../Pricing Library/Models')

% options' data
Prices_table = readtable('Example_Data.xlsx', 'Range', 'E1:G23');

mkt_prices = Prices_table.Price;
strikes    = Prices_table.Strike;
maturities = Prices_table.TimeToMaturity;


% plot ──────────────────────────────────────

red   = [204/255, 0/255, 0/255];
blue  = [21/255, 110/255, 235/255];
green = [40/255, 161/255, 0/255];

figure
set(gcf, 'Color', 'w', 'Name', 'Market Prices', 'NumberTitle', 'off')
hold on
grid on

plot(strikes, mkt_prices, ...
     'Color', red, ...
     'Linewidth', 1.2, ...
     'Marker', '.', 'Markersize', 12, 'MarkerEdgeColor', red)

xlabel('Strike : K')
ylabel('Prices')
title('Options'' prices')

%{
% put-call parity
mkt_prices(1:13) = mkt_prices(1:13) + S0 - strikes(1:13) .* exp( - r * maturities(1:13) );


% plot ──────────────────────────────────────

figure
set(gcf, 'Color', 'w', 'Name', 'Market Prices', 'NumberTitle', 'off')
hold on
grid on

plot(strikes(1:13), mkt_prices(1:13), ...
     'Color', blue, ...
     'Linewidth', 1.2, ...
     'Marker', '.', 'Markersize', 12, 'MarkerEdgeColor', blue)
plot(strikes(14:end), mkt_prices(14:end), ...
     'Color', green, ...
     'Linewidth', 1.2, ...
     'Marker', '.', 'Markersize', 12, 'MarkerEdgeColor', green);

xlabel('Strike : K')
ylabel('Prices')
title('Options'' prices')
%}


%% Calibration

% minimization parameters
LB = [eps;  -3; eps; eps];
x0 = [0.2; 0.5; 0.5; 0.2]; % sigma, theta, k, sigma bar
UB = [  1;   3;   3;   1];

% specify options for lsqnonlin
options = optimset('Display', 'off', ...
                   'PlotFcns', 'my_optimplotfval');
%{
options = optimset('Display', 'off', ...
                   'PlotFcns', 'my_optimplotfval'), ..
                   'TolFun', 1e-6, ...
                   'MaxIter', 1e6, ...
                   'MaxFunEvals', 1e6);
%}

% ┌----- PRICES -----------------------------------------------------------------------------------------------------------┐
%{
% parameters calibration
PARAMS = ...
    lsqnonlin(@(p) abs(arrayfun(@(i) callPrice_FFT_ExtVG(S0, maturities(i), r, strikes(i), p) ...
                                      - mkt_prices(i), [1:length(strikes)]')), x0, LB, UB, options);
%}
% └------------------------------------------------------------------------------------------------------------------------┘

% ┌----- VOLATILITY -------------------------------------------------------------------------------------------------------┐

% parameters calibration
PARAMS = ...
    lsqnonlin(@(p) abs(arrayfun(@(i) blsimpv(S0, strikes(i), r, maturities(i), callPrice_FFT_ExtVG(S0, maturities(i), r, strikes(i), p)) ...
                                      - blsimpv(S0, strikes(i), r, maturities(i), mkt_prices(i)), [1:length(strikes)]')), x0, LB, UB, options);

% └------------------------------------------------------------------------------------------------------------------------┘

% compute the calibrated model prices
prices = arrayfun(@(i) callPrice_FFT_ExtVG(S0, maturities(i), r, strikes(i), PARAMS), [1:length(strikes)]');

% compute implied volatilities
mkt_vol = blsimpv(S0, strikes, r, maturities, mkt_prices);
model_vol = blsimpv(S0, strikes, r, maturities, prices);

% calibration errors
ERROR_PRICES = sum(abs(prices - mkt_prices));
ERROR_VOL = sum(abs(model_vol - mkt_vol));


% ┌----- PLOTS ------------------------------------------------------------------------------------------------------------┐

% plot market prices VS model prices
figure
set(gcf, 'Color', 'w', 'Name', 'Market Prices VS Model Prices', 'NumberTitle', 'off')

plot(strikes, mkt_prices, 'sm', 'Markersize', 7, 'Linewidth', 1);
hold on
plot(strikes, prices, '+b', 'Markersize', 7, 'Linewidth', 1);

% display absolute errors on graph
text(strikes - 3, min(mkt_prices, prices) - 1.5, num2str(abs(mkt_prices - prices)), 'FontSize', 6); % error for each point
t = text(min(strikes), mkt_prices(end), ['\bf Total Error : ', num2str(ERROR_PRICES)]);             % total error
t.FontSize = 13;

legend('Market prices', 'Calibrated Extended VG prices')
xlabel('Strike : K')
ylabel('Prices')
title('Market Prices VS Model Prices')
grid on


% plot market volatility VS model volatility
figure
set(gcf, 'Color', 'w', 'Name', 'Market Implied Volatility VS Model Implied Volatility', 'NumberTitle', 'off')

plot(strikes, mkt_vol, 'sm', 'Markersize', 7, 'Linewidth', 1);
hold on
plot(strikes, model_vol, '+b', 'Markersize', 7, 'Linewidth', 1);

% display absolute errors on graph
text(strikes - 3, min(mkt_vol, model_vol) - 0.005, num2str(abs(mkt_vol - model_vol)), 'FontSize', 6); % error for each point
t = text(min(strikes), mkt_vol(end), ['\bf Total Error : ', num2str(ERROR_VOL)]);                     % total error
t.FontSize = 13;

legend('Market Implied Volatility', 'Extended VG Implied Volatility')
xlabel('Strike : K')
ylabel('Implied Volatility')
title('Market Implied Volatility VS Model Implied Volatility')
grid on

% └------------------------------------------------------------------------------------------------------------------------┘


end

