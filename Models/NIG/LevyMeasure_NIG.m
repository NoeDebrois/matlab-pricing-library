function LEVY_MEASURE = LevyMeasure_NIG(y, PARAMS)

% Levy measure
% Model: NIG

%% Parameters

% model parameters
SIGMA = PARAMS(1);
THETA = PARAMS(2);
K_NIG = PARAMS(3);


%% Levy measure

A = THETA / SIGMA^2;
B = sqrt( THETA^2 + SIGMA^2./K_NIG )/SIGMA^2;
C = sqrt( THETA^2 + SIGMA^2/K_NIG )/( pi * SIGMA * sqrt(K_NIG) );

LEVY_MEASURE = C .* exp( A .* y ) .* besselk(1, B .* abs(y) ) ./ ( abs(y) );


end

