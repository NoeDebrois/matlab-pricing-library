function S = simulate_ExtVG(S0, T, r, N_SIM, N, PARAMS)

% Simulate Extended Variance Gamma prices with Monte Carlo technique starting from Gamma subordinators
% S(t) = S0 * exp( rt + X(t) )
% X(t) is an Extended Variance Gamma process in the Risk Neutral measure

%% Parameters

% model parameters
SIGMA     = PARAMS(1);
THETA     = PARAMS(2);
K_VG      = PARAMS(3);
SIGMA_BAR = PARAMS(4);

% delta time
dt = T/N;


%% Computations

% > gamma simulation

a = dt/K_VG;
b = 1;

% 1. icdf simulation
%G = K_VG * icdf('Gamma', rand(N_SIM, N), a, b);

% 2. gamrnd simulation
G = K_VG * gamrnd(a, b, N_SIM, N);

% 3. Johnk's or Best's simulation
%{
if a <= 1

    % initialize
    X = zeros(N_SIM, N);
    Y = zeros(N_SIM, N);
    ONES = ones(N_SIM, N);
    
    % compute
    while any( reshape( X + Y <= ONES , N_SIM * N, 1) )
        U1 = unifrnd(0, 1, N_SIM, N);
        U2 = unifrnd(0, 1, N_SIM, N);
        indexes = ( X + Y <= ONES );
        X(indexes) = U1(indexes).^(1/a);
        Y(indexes) = U2(indexes).^(1/(1 - a));
    end
    E = icdf('Exponential', rand(N_SIM, N), a);

    % return
    G = K_VG * X .* E ./ ( X + Y );

elseif a > 1

    % initialize
    b = a - 1;
    c = 3 * a - 3/4;
    U1 = unifrnd(0, 1, N_SIM, N);
    U2 = unifrnd(0, 1, N_SIM, N);
    W = U1 .* (1 - U1);
    Y = sqrt(c./W) .* (U2 - 1/2);
    X = b + Y;
    Z = 64 .* W.^3 .* U2.^3;

    % compute
    while any( reshape( log(Z) <= 2 .* (b .* log(X./b) - Y) , N_SIM * N, 1) )
        U1 = unifrnd(0, 1, N_SIM, N);
        U2 = unifrnd(0, 1, N_SIM, N);
        indexes = ( log(Z) <= 2 .* (b .* log(X./b) - Y) );
        W(indexes) = U1(indexes) .* (1 - U1(indexes));
        Y(indexes) = sqrt(c./W(indexes)) .* (U2(indexes) - 1/2);
        X(indexes) = b + Y(indexes);
        if any( reshape(X > zeros(N_SIM, N), N_SIM * N, 1) )
            Z(indexes) = 64 .* W(indexes).^3 .* U2(indexes).^3;
        else
            Z(indexes) = exp( 2 .* (b .* log(X(indexes)./b) - Y) ) - 1000;
        end
    end

    % return
    G = K_VG * X;

end
%}

% > characteristic exponent
charexp = @(u) - SIGMA_BAR^2/2 * u.^2 - 1 / K_VG * log( 1 + 0.5 * u.^2 * SIGMA.^2 * K_VG - 1i * THETA * K_VG * u );
drift = r - charexp(-1i); % in order to be risk-neutral


% > compute logreturns

% initialization
X = zeros(N_SIM, N+1);

for i = 1:N
    
    X(:, i+1) = X(:, i) + drift * dt + SIGMA_BAR * sqrt(dt) * randn(N_SIM, 1) + THETA * G(:, i) + SIGMA * sqrt(G(:, i)) .* randn(N_SIM, 1);

end

S = S0 * exp( X );


end

