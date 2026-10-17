function [Phi, C] = GramSchmidt(S, dt)
    % Inputs:
    % S  : M x L matrix containing M signals, each with L time samples.
    % dt : The time step size used for numerical integration.
    %
    % Outputs:
    % Phi: N x L matrix containing the N orthonormal basis functions.
    % C  : M x N matrix containing the constellation coefficients.

    [M, L] = size(S);
    Phi = []; % Initialize an empty array to hold our basis functions
    
    % Step 1: Find the orthonormal basis functions
    for i = 1:M
        g = S(i, :); % Grab the current signal
        
        % Subtract the projection of the current signal onto all existing basis functions
        if ~isempty(Phi)
            for j = 1:size(Phi, 1)
                % Projection integral: integral( S(i,t) * Phi(j,t) dt )
                projection = sum(S(i, :) .* Phi(j, :)) * dt;
                g = g - projection * Phi(j, :);
            end
        end
        
        % Calculate the energy of the remaining orthogonal component
        Eg = sum(g.^2) * dt;
        
        % If the energy is greater than a tiny threshold, it's a valid new basis function
        % (We use 1e-6 instead of 0 to account for MATLAB's floating-point inaccuracies)
        if Eg > 1e-6 
            phi_new = g / sqrt(Eg); % Normalize by dividing by sqrt of Energy
            Phi = [Phi; phi_new];   % Append to our basis matrix
        end
    end
    
    % Step 2: Calculate the Constellation Coefficients (C matrix)
    N = size(Phi, 1);
    C = zeros(M, N);
    
    for i = 1:M
        for j = 1:N
            % Find how much of basis 'j' is present in signal 'i'
            C(i, j) = sum(S(i, :) .* Phi(j, :)) * dt;
        end
    end
end