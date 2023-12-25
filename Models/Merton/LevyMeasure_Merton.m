function LEVY_MEASURE = LevyMeasure_Merton(y, PARAMS)

% Levy measure
% Model: Merton

%% Parameters

% model parameters
SIGMA  = PARAMS(1);
LAMBDA = PARAMS(2);
MU     = PARAMS(3);
DELTA  = PARAMS(4);


%% Levy measure

LEVY_MEASURE = LAMBDA * exp( - ( y - MU ).^2/( 2 * DELTA^2 ) )/sqrt( 2 * pi * DELTA^2 );


end

