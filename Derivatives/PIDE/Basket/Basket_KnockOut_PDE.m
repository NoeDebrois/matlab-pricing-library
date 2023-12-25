

clc
clear all
close all

set(0, 'DefaultFigureWindowStyle', 'docked');

% Pricing a 2D basket option (Knock & Out Call)

%% Parameters

% option parameters
T    = 1;
S_10 = 100;
S_20 = 70;
K    = 190;

% barrier parameters
L1 = 80;
U1 = 120;
L2 = 60;
U2 = 110;

% model parameters
r = 0.01;
sigma1 = 0.4;
sigma2 = 0.35;
rho = - 0.3;

% discretization parameters
MT = 50; % points in time
M1 = 30; % points in x1
M2 = 40; % points in x2


%% Grid

% space grid
x1_min = log(L1); x1_max = log(U1);
x2_min = log(L2); x2_max = log(U2);
x1 = linspace(x1_min, x1_max, M1); dx1 = x1(2) - x1(1);
x2 = linspace(x2_min, x2_max, M2); dx2 = x2(2) - x2(1);
[X1, X2] = meshgrid(x1, x2);
X = [X1(:) X2(:)];

% time grid
dt = T/MT;

% total number of points
numpoints = M1 * M2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative scheme :
%
%            i+1
%  i-M2       i      i+M2
%            i-1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% boundary indexes
W = [1:M2];
S = [1:M2:M1*M2];
N = [M2:M2:M1*M2];
E = [(M1-1)*M2+1:M1*M2];
boundary_idx = unique([ W, S, N, E ]);


%% Matrix

M = sparse(numpoints, numpoints);

for i = 1:numpoints

    if min( abs( i - boundary_idx ) ) == 0

        M(i, i) = 1;

    else

        % derivative wrt time and -rV term
        M(i, i) = - 1/dt - r;
        
        % first order derivative wrt x1
        coeff = r - sigma1^2/2;
        M(i, i + M2) = coeff/( 2 * dx1 );
        M(i, i - M2) = - coeff/( 2 * dx1 );
        
        % first order derivative wrt x2
        coeff = r - sigma2^2/2;
        M(i, i + 1) = coeff/( 2 * dx2 );
        M(i, i - 1) = -coeff/( 2 * dx2 );
        
        % second order derivative wrt x1
        coeff = sigma1^2/2;
        M(i, i + M2) = M(i, i + M2) + coeff/dx1^2;
        % Note: since we are adding a coefficient to the one already
        % present in M(i, i + M2) we do it as M(i, i + M2) = M(i, i + M2) + ...
        M(i, i) = M(i, i) - 2 * coeff/dx1^2;
        M(i, i - M2) = M(i, i - M2) + coeff/dx1^2;
        
        % second order derivative wrt x2
        coeff = sigma2^2/2;
        M(i, i + 1) = M(i, i + 1) + coeff/dx2^2;
        M(i, i) = M(i, i) - 2 * coeff/dx2^2;
        M(i, i - 1) = M(i, i - 1) + coeff/dx2^2;
        
        % mixed second order derivative
        coeff=rho*sigma1*sigma2;
        M(i, i + M2 + 1) = coeff/( 4 * dx1 * dx2 );
        M(i, i - M2 + 1) = -coeff/( 4 * dx1 * dx2 );
        M(i, i + M2 - 1) = - coeff/( 4 * dx1 * dx2 );
        M(i, i - M2 - 1) = coeff/( 4 * dx1 * dx2 );

    end

end


%% Implicit Euler

% payoff at expiry (call)
V = max( exp( X(:, 1) ) + exp( X(:, 1) ) - K, 0);

% boundary conditions (knock-out call)
V(boundary_idx) = 0;

for j = MT:-1:1

    V = M \ ( - V / dt );

end


%% Output

figure
surf(exp(X1), exp(X2), reshape(V, size(X1)))

% price
price = griddata(exp(X1), exp(X2), reshape(V, size(X1)), S_10, S_20);

fprintf('\nKnock-Out basket option\n')
fprintf('Price:  %.5f\n\n', price)

