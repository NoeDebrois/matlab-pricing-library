function LEVY_MEASURE = LevyMeasure_VG(y, PARAMS)

% Levy measure
% Model: VG

%% Parameters

% model parameters
SIGMA     = PARAMS(1);
THETA     = PARAMS(2);
K_VG      = PARAMS(3);


%% Levy measure

A = THETA / SIGMA^2;
B = sqrt( THETA^2 + 2 * SIGMA^2/K_VG )/SIGMA^2;

LEVY_MEASURE = exp( A * y - B * abs(y) ) ./ ( K_VG * abs(y) );


end

