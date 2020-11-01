function D = homogenise2DSubproblem( subproblem, boundary_type, skin_map ) 
% This function takes an input subproblem (specified in the same format as
% a full problem object) and solves it for two different imposed gradients,
% thus providing the requisite information for determining the effective
% conductivity tensor. A skin map can also be supplied, which specifies
% which sites should or should not be used in calculating the average
% fluxes (the full region specified by 'subproblem' will still be the
% region over which the closure problem is solved).

% Specifies at what point tensor elements may be considered essentially
% zero (values smaller than this in magnitude will be set to zero)
zero_tol = 1e-12;

% Specify whether to display warnings or not when diffusion tensors are not
% SPD and are forcibly corrected
show_warnings = false;

% Default to no skin if no variable skin_map was provided
if nargin < 2
    skin_map = false(size(subproblem.occ_map));   % Just use the occupancy map as an easy way to read out the size of the input subproblem
end
    

% Solve the closure problem for imposed gradients in the 'x' and 'y' directions
switch boundary_type
    
    case {'linear','confined'}
        [Jx_bar_x, Jy_bar_x] = solveClosureProblem( subproblem, 'x', boundary_type, skin_map);
        [Jx_bar_y, Jy_bar_y] = solveClosureProblem( subproblem, 'y', boundary_type, skin_map);
        
    case 'periodic'
        [Jx_bar_x, Jy_bar_x] = solveClosureProblemPeriodic( subproblem, 'x', skin_map );
        [Jx_bar_y, Jy_bar_y] = solveClosureProblemPeriodic( subproblem, 'y', skin_map );
        
    case 'hybrid'
        [Jx_bar_x, ~] = solveClosureProblem( subproblem, 'x', 'confined', skin_map );
        [~, Jy_bar_y] = solveClosureProblem( subproblem, 'y', 'confined', skin_map );
        [~, Jx_bar_y] = solveClosureProblem( subproblem, 'x', 'linear', skin_map );
        [Jy_bar_x, ~] = solveClosureProblem( subproblem, 'y', 'linear', skin_map );
        
    otherwise
        error('These type of homogenisation boundary conditions have not been implemented. Specify either ''linear'', ''confined'', ''periodic'', or ''hybrid''.');
        
end

% Use these averaged fluxes as the conductivity tensor elements
D_xx = Jx_bar_x;
D_xy_x = Jy_bar_x;
D_xy_y = Jx_bar_y;
D_yy = Jy_bar_y;

% Forcibly symmeterise the conductivity tensor.
D_xy = ( D_xy_x + D_xy_y ) / 2;

% Convert the individual values to a matrix
D = [ [D_xx, D_xy]; [D_xy, D_yy] ];

% Zero out any values that are sufficiently close to zero to disregard
D( abs(D) < zero_tol ) = 0;

% Check if the diffusion tensor is symmetric positive definite. If not,
% zero out the negative eigenvalues and rebuild it
[v, lambda] = eig(D);
if any(diag(lambda) < 0)
    
    % Grab out eigenvalues
    lambda = diag(lambda);
    
    % Check for any being complex - error if so
    if any(~isreal(lambda))
        error('Matrix has complex eigenvalues!');
    end
    
    % Zero any negative eigenvalues
    lambda(lambda < 0) = 0;
    
    % Rebuild the matrix
    D = v * diag(lambda) * inv(v);
    
    % Output a warning in this circumstance too (if flag set)
    if show_warnings
        warning('Output diffusion tensor was forcibly made SPD!'); 
    end
end