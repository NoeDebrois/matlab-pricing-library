function [S, Sav] = simulate_Bates_AV(S0, T, r, N_SIM, N, PARAMS)

% Implementation of QE Scheme for Bates model with antithetic variables
% "Efficient Simulation of the Heston Stochastic Volatility Model", Andersen

%% Parameters

epsilon = PARAMS(1); % vol-of-vol
k       = PARAMS(2); % mean reversion speed
rho     = PARAMS(3); % correlation
theta   = PARAMS(4); % mean
V0      = PARAMS(5); % starting value of the variance

% jumps parameters
mu      = PARAMS(6);
delta   = PARAMS(7);
lambda  = PARAMS(8);
% Note: the starting point of the variance is not known, so we have to calibrate it

dt = T/N; 
V  = zeros(N_SIM, N+1);
S  = zeros(N_SIM, N+1);
Sav = zeros(N_SIM, N+1);

% Random changes for volatility are sampled from one of two distributions,
% depending on the ratio Psi = s^2/m^2, where m & s are mean and variance
% of next volatility value, conditioned to current one.
%  Scheme 1: exponential approximation - selected when Psi > Psi_cutoff
%  Scheme 2: quadratic aapproximation  - selected when 0 < Psi < Psi_cutoff

% Psi_cutoff
% Note: choose Psi_cutoff = 1.5 as suggested in article
Psi_cutoff = 1.5;


%% Discetize V

V(:, 1) = V0;

for i = 1:N

    % STEP 1 & 2
    % -> calculate m, s, Psi

    m  = theta + ( V(:, i) - theta ) * exp( - k * dt );
    m2 = m .^ 2;
    s2 = V(:, i) * epsilon ^2 * exp( - k * dt ) * ( 1 - exp( - k * dt ) ) / k ...
          + theta * epsilon ^2 * ( 1 - exp( - k * dt ) ) ^2 / ( 2 * k );
    s  = sqrt( s2 );
    
    Psi = ( s2 ) ./ ( m2 );
    

    % STEP 3,4,5
    % -> depending on Psi, use exp or quad scheme to calculate next V

    % 1. Exponential approximation - used for Psi > Psi_cutoff
    % The PDF of V(t+dt) is p * delta(0) + (1 - p) * (1 - exp( - beta x )
    % thus a probability mass in 0 and an exponential tail after that
    
    % take the simulations for which Psi > Psi_cutoff
    index = find( Psi > Psi_cutoff );

    p_exp = ( Psi(index) - 1 )./( Psi(index) + 1 ); % probability mass in 0                - Eq 29
    beta_exp = ( 1 - p_exp )./m(index);             % exponent of exponential density tail - Eq 30
    
    % gets x from inverse CDF applied to uniform U
    U = rand(size(index));
    V(index, i+1) = ( log( ( 1 - p_exp ) ./ ( 1 - U ) ) ./ beta_exp ) .* ( U > p_exp );

    % 2. Quadratic approx - used for 0 < Psi < Psi_cutoff
    % V(t+dt) = a( b + Zv )^2, Zv ~ N(0, 1)
    
    % take the simulations for which 0 < Psi < Psi_cutoff
    index = find( Psi <= Psi_cutoff );

    invPsi  = 1 ./ Psi(index);
    b2_quad = 2 * invPsi - 1 + sqrt( 2 * invPsi ) .* sqrt( 2 * invPsi - 1 );    %  - Eq 27
    a_quad  = m(index)./( 1 + b2_quad );                                        %  - Eq 28
    V(index, i+1) = a_quad .* ( sqrt( b2_quad ) + randn( size( index ) ) ) .^2; %  - Eq.23
    
end


%% Discetize X

S(:, 1)   = S0;
Sav(:, 1) = S0;
gamma1 = 0.5; % => central discretization scheme
gamma2 = 0.5; % => central discretization scheme

% fix drift to obtain risk-neutrality
CharExp = @(v) lambda * ( exp( - delta^2 * v.^2/2 + 1i * mu * v ) - 1 ); % jumps characteristic exponent without drift
drift_rn = - CharExp(-1i);                                               % drift chosen under the risk neutral measure

k0 = ( r + drift_rn ) * dt - rho * k * theta * dt / epsilon;
k1 = gamma1 * dt * ( k * rho / epsilon - 0.5 ) - rho / epsilon;
k2 = gamma2 * dt * ( k * rho / epsilon - 0.5 ) + rho / epsilon;
k3 = gamma1 * dt * ( 1 - rho^2 );
k4 = gamma2 * dt * ( 1 - rho^2 );

% Gaussian random variables
Z = randn(N_SIM, N);

% number of jumps
NT = icdf('Poisson', rand(N_SIM, 1), lambda * T);


for j = 1:N_SIM
    
    JumpTimes = sort( T * rand(NT(j), 1) );
    
    for i = 1:N
        
        % add diffusion component
        S(j, i+1)   = exp( log( S(j, i) ) ...
                      + k0                ...
                      + k1 * V(j, i)      ...
                      + k2 * V(j, i+1)    ...
                      + sqrt( k3 * V(j, i) + k4 * V(j, i+1) ) .* Z(j, i) );
        Sav(j, i+1) = exp( log( Sav(j, i) ) ...
                      + k0                ...
                      + k1 * V(j, i)      ...
                      + k2 * V(j, i+1)    ...
                      - sqrt( k3 * V(j, i) + k4 * V(j, i+1) ) .* Z(j, i) );

        % add jump part if exists jump in ((i-1)*dt, i*dt)
        for l = 1:NT(j)
            
            if ( ( JumpTimes(l) > (i - 1) * dt ) && ( JumpTimes(l) <= i * dt ) )
                
                Z_jump = randn;
                Y = mu + delta * Z_jump;
                Yav = mu - delta * Z_jump;
             
                S(j, i+1)   = S(j, i+1)   * exp( Y );
                Sav(j, i+1) = Sav(j, i+1) * exp( Yav );
            
            end
            
        end
         
    end
    
end


end

