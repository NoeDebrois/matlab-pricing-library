function [S, Sav] = simulate_Kou_AV(S0, T, r, N_SIM, N, PARAMS)

% Simulate Kou prices with Monte Carlo antithetic variables technique starting from the logreturns X(t)
% S(t) = S0 * exp( rt + X(t) )
% X(t) is a Kou process in the Risk Neutral measure

%% Parameters

% model parameters
SIGMA        = PARAMS(1);
LAMBDA       = PARAMS(2);
LAMBDA_MINUS = PARAMS(3);
LAMBDA_PLUS  = PARAMS(4);
P            = PARAMS(5);

% delta time
dt = T/N;


%% Computations

% initialize X and Xav
X = zeros(N_SIM, N+1); % X = rt + X(t)
Xav = X;

% number of jumps
NT = icdf('Poisson', rand(N_SIM, 1), LAMBDA * T);

% characteristic exponent
charexp = @(u) - SIGMA^2 * u.^2/2 + 1i * u * LAMBDA .* ( P ./ ( LAMBDA_PLUS - 1i * u ) - ( 1 - P ) ./ ( LAMBDA_MINUS + 1i * u ) );
drift = r - charexp( - 1i ); % in order to be risk-neutral

% sample standard normal
Z = randn(N_SIM, N);


for j = 1:N_SIM

    JumpTimes = sort( T * rand(NT(j), 1) );

    for i = 1:N

        % add diffusion component
        X(j, i+1)   = X(j, i) + drift * dt + SIGMA * sqrt(dt) * Z(j, i);
        Xav(j, i+1) = X(j, i) + drift * dt - SIGMA * sqrt(dt) * Z(j, i);

        % add jump part -> ( (i-1)dt, idt )
        for l = 1:NT(j)

            if JumpTimes(l) > (i-1) * dt && JumpTimes(l) <= i * dt

                sim_p = rand;

                if sim_p < P % positive jump
                    Y = icdf('exp', rand, 1/LAMBDA_PLUS);
                else % negative jump
                    Y = - icdf('exp', rand, 1/LAMBDA_MINUS);
                end

                X(j, i+1)   = X(j, i+1)   + Y;
                Xav(j, i+1) = Xav(j, i+1) + Y;
                % Note: the antithetic variables are added only to the brownian motion, not Y

            end
            
        end

    end

end

S = S0 * exp( X );
Sav = S0 * exp( Xav );


end

