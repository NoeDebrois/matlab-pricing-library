function [S, Sav] = simulate_Merton_AV(S0, T, r, N_SIM, N, PARAMS)

% Simulate Merton prices with Monte Carlo antithetic variables technique starting from the logreturns X(t)
% S(t) = S0 * exp( rt + X(t) )
% X(t) is a Merton process in the Risk Neutral measure

%% Parameters

% model parameters
SIGMA  = PARAMS(1);
LAMBDA = PARAMS(2);
MU     = PARAMS(3);
DELTA  = PARAMS(4);

% delta time
dt = T/N;


%% Computations

% initialize X
X = zeros(N_SIM, N+1); % X = rt + X(t)
Xav = X;

% number of jumps
NT = icdf('Poisson', rand(N_SIM, 1), LAMBDA * T);

% characteristic exponent
charexp = @(u) - 0.5 * SIGMA^2 * u.^2 + LAMBDA * ( exp(- 0.5 * DELTA^2 * u.^2 + 1i * MU * u) - 1 );
drift = r - charexp(-1i); % in order to be risk-neutral

% sample standard normal
Z = randn(N_SIM, N);


for j = 1:N_SIM
    
    JumpTimes = sort( T * rand(NT(j), 1) );
    
    for i = 1:N

        % add diffusion component
        X(j, i+1)   = X(j, i)   + drift * dt + SIGMA * sqrt(dt) * Z(j, i);
        Xav(j, i+1) = Xav(j, i) + drift * dt - SIGMA * sqrt(dt) * Z(j, i);
        
        % add jump part if exists jump in ((i-1)dt, idt)
        for l = 1:NT(j)
            
            if JumpTimes(l) > (i - 1) * dt && JumpTimes(l) <= i * dt
            
                g = randn;
                Y   = MU + DELTA * g;
                Yav = MU - DELTA * g;
                
                X(j, i+1)   = X(j, i+1)   + Y;
                Xav(j, i+1) = Xav(j, i+1) + Yav;
            
            end
            
        end
         
    end
    
end

S   = S0 * exp( X );
Sav = S0 * exp( Xav );


end

