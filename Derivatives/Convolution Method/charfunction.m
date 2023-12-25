function F = charfunction(u, parameters, flag)

% flag = 0 => characteristic function for the backward problem
% flag = 1 => characteristic function for the forward problem

if nargin == 2
    flag = 0; % characteristic function for the backward problem by default
end

% risk-neutral characteristic exponent
meancorrection = (parameters.rf - parameters.q) * parameters.dt - log( charfunction0( -1i, parameters) );
         
% compute characteristic function
F = exp(1i * meancorrection * u) .* charfunction0(u, parameters);

if flag == 0
    F = conj(F);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = charfunction0(u, parameters)

dt = parameters.dt;

switch parameters.distr

    case 1 % Normal
        
        % parameters
        m = parameters.m;
        s = parameters.s;
	
	    % rearrange parameters (time rescaling)
        m = m * dt;
        s = s * sqrt(dt);

        % compute characteristic function
        F = exp(1i * u * m - 0.5 * (s * u).^2);

    case 2 % Normal inverse Gaussian (NIG)
        
        % parameters
        alpha = parameters.alpha;
        beta = parameters.beta;
        delta = parameters.delta;

        % rearrange parameters (time rescaling)
        alpha = alpha;
        beta = beta;
        delta = delta * dt;
        
        % compute characteristic function
        F = exp( - delta * ( sqrt( alpha^2 - ( beta + 1i * u).^2 ) - sqrt( alpha^2 - beta^2 ) ) );
   
end


end
