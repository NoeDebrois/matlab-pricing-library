function LEVY_MEASURE = LevyMeasure_Kou(y, PARAMS)

% Levy measure
% Model: Kou

%% Parameters

% model parameters
SIGMA        = PARAMS(1);
LAMBDA       = PARAMS(2);
LAMBDA_MINUS = PARAMS(3);
LAMBDA_PLUS  = PARAMS(4);
P            = PARAMS(5);


%% Levy measure

LEVY_MEASURE = P .* LAMBDA .* LAMBDA_PLUS .* exp( - LAMBDA_PLUS .* y ) .* ( y > 0 ) ...
                + ( 1 - P ) .* LAMBDA .* LAMBDA_MINUS .* exp( - LAMBDA_MINUS .* abs( y ) ) .* ( y < 0 );


end

