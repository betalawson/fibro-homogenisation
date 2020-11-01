function homog_problem = homogeniseFull2DProblemRetainGrid(problem, Ax, Ay, boundary_type)
% This function takes an input problem, then uses block homogenisation to
% reduce its number of mesh points. Ax and Ay are the number of smaller
% elements (in each direction) to absorb into the larger elements. In the
% current code, Ax and Ay must divide without remainder into the original
% numbers of elements in each dimension, because the solver code at current
% handles only regular grids

% Specify if 'skin' is to be used - borrowing information from surrounding
% blocks in order to improve determination of effective diffusivity and
% also reduce emphasis on selection of boundary conditions. Skin will be
% applied using half of the size of a single averaging volume on each edge
use_skin = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read out some basic details of the problem
occ_map = problem.occ_map;

% Read out the number of microscopic volumes
[Ny, Nx] = size(occ_map);

% Determine the number of averaging volumes that will be constructed
% (according to specified dimensions)
NAx = ceil(Nx / Ax);
NAy = ceil(Ny / Ay);

% Initialise the homogenised problem to the original problem (because in
% this case, not much will change)

if use_skin
    
    % Determine the amount of skin to use. For now, half of the block size
    % is chosen, in accordance with the literature
    Sx = ceil(Ax/2);
    Sy = ceil(Ay/2);
    
    % Create a base skin map, which is the skin map to use if all border
    % directions include skin
    skin_map = [true(Sy,Ax + 2*Sx); [true(Ay, Sx), false(Ay,Ax), true(Ay, Sx)]; true(Sy,Ax + 2*Sx)];
    
end


% Create a new set of macroscopic volumes, of the specified dimension
Aocc = false(Ny, Nx); AVfrac = zeros(Ny, Nx);
AD_xx = zeros(Ny, Nx);   AD_xy = zeros(Ny, Nx);   AD_yy = zeros(Ny, Nx);

for i = 1:NAx
    for j = 1:NAy
        
        % Define the regions in the microscopic matrix that the current
        % averaging volume occupies
        iV = (i-1) * Ax + (1:Ax);
        jV = (j-1) * Ay + (1:Ay);
        
        % Check if homogenisation is actually required here (not required
        % if a region is completely blocked). First, read out the occupancy
        % (discluding any skin)
        occ_here = occ_map(jV,iV);
        
        % Use this occupancy map to calculate the volume fraction
        Vfrac_here = sum(~occ_here(:)) / numel(occ_here);
        AVfrac(jV,iV) = Vfrac_here;
        
        if Vfrac_here == 0    % If all sites are occupied in this averaging volume, set diffusivity to zero
            
            % This is a completely occupied site
            Aocc(jV,iV) = true;
            % All diffusion here is zero
            AD_xx(jV,iV) = 0;
            AD_xy(jV,iV) = 0;
            AD_yy(jV,iV) = 0;
            
        else   % Otherwise, we need to homogenise
            
            % Add skin if requested
            if use_skin
                
                % Assume that skin will be included in all directions, but
                % then mark sites that cannot be included to be removed
                usable_sites = true( Ay + 2*Sy, Ax + 2*Sx );
                
                if i == 1   % Cannot use skin on left (or corners on left)
                    usable_sites(:, 1:Sx) = false;
                end
                if i == NAx  % Cannot use skin on right (or corners on right)
                    usable_sites(:, Ax+Sx+(1:Sx) ) = false;
                end
                if j == 1   % Cannot use skin on bottom (or corners on bottom)
                    usable_sites(1:Sy, :) = false;
                end
                if j == NAy  % Cannot use skin on top (or corners on top)
                    usable_sites(Ay+Sy+(1:Sy), :) = false;
                end
                
                % Use the portions of the default skin map that are marked
                % as used sites
                skin_map_here = skin_map(usable_sites');
                
                % To obtain the sites that will be used, pretend a full set
                % will be marked as true, but only actually mark as true
                % those sites which were also marked usable. This will
                % prevent any elements outside of the bounds of the matrix
                % from being attempted to be set to true. In order to
                % achieve this, first a list of indices for the sites
                % to fill in must be created
                [index_i, index_j] = meshgrid( [iV(1)-Sx-1 + (1:Sx), iV, iV(end) + (1:Sx)], [jV(1)-Sy-1 + (1:Sy), jV, jV(end) + (1:Sy)] );
                index_j = index_j(usable_sites);
                index_i = index_i(usable_sites);
                
                % Grab out the range of sites to actually use
                uiV = unique(index_i);
                ujV = unique(index_j);
                
            else
                
                % If not using skin, just set the skin map to empty so it's
                % not used by the homogenisation code (jV and iV work
                % unchanged)
                skin_map_here = [];
                uiV = iV;
                ujV = jV;
                
            end
            
            % Set up the sub problem
            
            
            subproblem.occ_map = occ_map(ujV,uiV);
            subproblem.D_tensor.D_xx = problem.D_tensor.D_xx(ujV,uiV);
            subproblem.D_tensor.D_xy = problem.D_tensor.D_xy(ujV,uiV);
            subproblem.D_tensor.D_yy = problem.D_tensor.D_yy(ujV,uiV);
            subproblem.Vfrac = problem.Vfrac(ujV,uiV);
            subproblem.grid = problem.grid;
            
            % Run the homogenisation code to calculate the diffusion tensor
            % for this block
            D_here = homogenise2DSubproblem( subproblem, boundary_type, skin_map_here );
            
            % Store its elements in the new homogenised problem structure
            AD_xx(jV,iV) = D_here(1,1);
            AD_xy(jV,iV) = D_here(1,2);
            AD_yy(jV,iV) = D_here(2,2);
            
            % This site is unoccupied (completely) by blockage
            Aocc(jV,iV) = false;
            
        end
        
    end
end


% The homogenised problem is just the base problem, but with the new
% diffusion tensors, volume fractions and occupancy map
homog_problem = problem;
homog_problem.occ_map = Aocc;
homog_problem.D_tensor.D_xx = AD_xx;
homog_problem.D_tensor.D_xy = AD_xy;
homog_problem.D_tensor.D_yy = AD_yy;
homog_problem.Vfrac = AVfrac;

end

