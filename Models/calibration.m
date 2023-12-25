
% bounds for parameters

- standard deviations (e.g. sigma, delta, eta)
Bounds : [eps, 1] or [0, 1] since a volatility of 100% is very high and a negative volatility is not admissible
Start  : around 0.2

- mean (e.g. mu, theta)
Bounds : [-5, 5] to have a symmetric interval
Start  : 0 to not give any advantage to having a positive or negative drift

- vol-of-vol (e.g. k)
Bounds : [eps, 5 * mean(maturities)] since k should be around the maturities of the options considered and it is positive
Start  : mean(maturities) since k should be around the maturities of the options considered

- jump intensities (e.g. lambda)
Bounds : [0, 200]
Start  : around 40

- parameter of exponential distributions (e.g. lambda plus, lambda minus)
Bounds : [0, 20] to have the biggest sensical interval
Start  : around 2 to reach sensical values

- probabilities (e.g. p)
Bounds : [0, 1] since are the natural bounds of the probability
Start  : 0.5 since starting in the middle does not give advantages to having a positive or a negative jump component

- mean reverting speed (e.g. xi)
Bounds : [eps, 3]
Start  : around 0.2


% how to decide between different starting points

- choose the one with the smallest error
- choose the one with the smallest maximum error
- look at the graph


% put-call parity

Data_put(:, 1) = Data_put(:, 1) + S0 - Data_put(:, 2) * exp( - r * maturity );
% call = put + S0 - K * exp( - r * T )


% how to understand if the data are from put or call

% European options with maturity close to 1 year (16 sept 2022)
Data_put = [1.83    95      236 % prices (USD), strikes (USD), volume
            2.01    97.5    100
            2.23   100      505
            2.45   102.5    142
            2.70   105      295
            2.99   107.5     79
            3.38   110      181
            3.70   112.5     16
            4.10   115      901
            4.50   117.5    454
            5.00   120      228
            5.60   122.5    185
            6.20   125     1253];

Data_call = [27.15   130    110 % prices (USD), strikes (USD), volume
             23.85   135    108
             20.75   140     88
             17.85	 145    388
             15.50	 150   1200
             13.23	 155    111
             11.40	 160    394
              9.72	 165     59
              8.21	 170    206
              7.00	 175     58
              5.88	 180    147];

mkt_prices = [Data_put(:, 1); Data_call(:, 1)]; % market prices
strikes    = [Data_put(:, 2); Data_call(:, 2)]; % strikes

figure
set(gcf, 'Color', 'w');
plot(strikes(1:13), mkt_prices(1:13), 's', 'Color', 'g', 'Markersize', 7, 'Linewidth', 1.5);
hold on
plot(strikes(14:end), mkt_prices(14:end), 'x', 'Color', 'r', 'Markersize', 7, 'Linewidth', 1.5);

legend('Put', 'Call')



