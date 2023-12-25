function [S, v] = CONV(S_0, K, Ndate, N, Barrier, param)

% Price a down & out call option using the convolution method

% xmin - xmax
b = 2.5;

% kernel
[x, ~, ~, H] = kernel(N, -b, b, param, 0);

% stock price
S = S_0 * exp(x);

% v(x, 0) = payoff function
v = max(S - K, 0) .* (S > Barrier);

% Fourier transform of the density
H = ifftshift(H); % swap the left and right halves of H


% CONVOLUTION ALGORITHM -
% Note: for us the ifft is the Fourier transform and fft is the inverse Fourier
% transform, becuase we take it with different signs with respect to MATLAB's ones
for j = 1:Ndate
    
    % Fourier transform of v
    FT_v = ifft( ifftshift(v) ); % swap the left and right halves of v
    
    % anti Fourier tranform of FT_v * H to find v at the next step
    v = real( fftshift( fft( FT_v .* H ) ) );
    
    % add the discount to find v_EU
    v_EU = v * exp( - param.rf * param.dt );
    
    % add the barrier
    v = v_EU .* (S > Barrier); % put equal to zero the values where S is smaller than the Barrier

end


% focus on the grid
index = find( ( S > 0.1 * S_0 ) .* ( S < 3 * S_0 ) );
S = S(index);
v = v(index);

% plot the values of the option for the different stock prices
figure
plot(S, v)


end

