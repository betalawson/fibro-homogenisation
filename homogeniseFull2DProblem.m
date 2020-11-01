function homog_problem = homogeniseFull2DProblem(problem, Ax, Ay, boundary_type)
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


if use_skin
    
    % Determine the amount of skin to use. For now, half of the block size
    % is chosen, in accordance with the literature
    Sx = ceil(Ax/2);
    Sy = ceil(Ay/2);
    
    % Create a base skin map, which is the skin map to use if all border
    % directions include skin
    skin_map = [true(Sy,Ax + 2*Sx); [true(Ay, Sx), false(Ay,Ax), true(Ay, Sx)]; true(Sy,Ax + 2*Sx)];
    
end


% Nodes will be placed at averaging volume boundaries (vertex-centred
% finite volume).  So first create node positions
[AnodeX, AnodeY] = meshgrid( linspace(0,problem.grid.Lx,NAx+1), linspace(0,problem.grid.Ly,NAy+1) );

% Create a new set of macroscopic volumes, of the specified dimension
Aocc = false(NAy, NAx); AVfrac = zeros(NAy, NAx);
AD_xx = zeros(NAy, NAx);   AD_xy = zeros(NAy, NAx);   AD_yy = zeros(NAy, NAx);
model_warning_flag = false;
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
        AVfrac(j,i) = sum(~occ_here(:)) / numel(occ_here);
        
        if AVfrac(j,i) == 0    % If all sites are occupied in this averaging volume, set diffusivity to zero
            
            % This is a completely occupied site
            Aocc(j,i) = true;
            % All diffusion here is zero
            AD_xx(j,i) = 0;
            AD_xy(j,i) = 0;
            AD_yy(j,i) = 0;
            
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
                    % as used sites - this also vectorises, so use the
                    % transpose of the usable sites matrix so that the
                    % correct order for vectors is maintained
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
                    iV = unique(index_i);
                    jV = unique(index_j);
                
            else
                
                % If not using skin, just set the skin map to empty so it's
                % not used by the homogenisation code (jV and iV work
                % unchanged)
                skin_map_here = [];
                
            end
            
            % Set up the sub problem.       
            subproblem.occ_map = occ_map(jV,iV);
            subproblem.D_tensor.D_xx = problem.D_tensor.D_xx(jV,iV);
            subproblem.D_tensor.D_xy = problem.D_tensor.D_xy(jV,iV);
            subproblem.D_tensor.D_yy = problem.D_tensor.D_yy(jV,iV);
            subproblem.Vfrac = problem.Vfrac(jV,iV);
            subproblem.grid = problem.grid;
            
            % Run the homogenisation code to calculate the diffusion tensor
            % for this block
            D_here = homogenise2DSubproblem( subproblem, boundary_type, skin_map_here );
            
            % Store its elements in the new homogenised problem structure
            AD_xx(j,i) = D_here(1,1);
            AD_xy(j,i) = D_here(1,2);    % Should be symmetric
            AD_yy(j,i) = D_here(2,2);
            
            % This site is not completely occupied by blockage
            Aocc(j,i) = false;
            
        end
        
    end
end


% Now perform a similar loop, but this time over nodes
Astim_sites1 = false(NAy+1, NAx+1);
Astim_sites2 = false(NAy+1, NAx+1);
Amodel_assignments = zeros(NAy+1, NAx+1);
for i = 1:NAx+1
    for j = 1:NAy+1
        
        % Read out micro nodes associated with this macro node
        inodesV = (i-1) * Ax + 1+( ( -floor(Ax/2)*(i>1) ): (floor(Ax/2) * (i <NAx+1)) );
        jnodesV = (j-1) * Ay + 1+( ( -floor(Ay/2)*(j>1) ): (floor(Ay/2) * (j <NAy+1)) );
        
        % Select which cell model to assign here by checking how many times
        % each model occurs in this averaging volume
        models_here = problem.model_assignments(jnodesV,inodesV);
        model_counts = zeros(length(problem.cell_models),1);
        for k = 1:length(problem.cell_models)
            model_counts(k) = sum(sum(models_here == k));
        end
        [max_count, max_model] = max(model_counts);
        Amodel_assignments(j,i) = max_model;
        if max_count < numel(models_here)
            model_warning_flag = true;
        end
        
        % To elect whether to assign the corner nodes as stimulus sites,
        % check if any nodes within the averaging volume's 'sub corners'
        % are stimulus sites
        if any(any( problem.stim_sites1(jnodesV,inodesV) ) )
            Astim_sites1(j,i) = true;
        end
        if any(any( problem.stim_sites2(jnodesV,inodesV) ) )
            Astim_sites2(j,i) = true;
        end
        
    end
end


if model_warning_flag
    fprintf('WARNING: One or more averaging volumes was composed of nodes with different cell models. The most frequently occurring cell model has been used. \n');
end


% Now create a new homogenised problem
homog_problem.occ_map = Aocc;
homog_problem.D_tensor.D_xx = AD_xx;
homog_problem.D_tensor.D_xy = AD_xy;
homog_problem.D_tensor.D_yy = AD_yy;
homog_problem.Vfrac = AVfrac;
homog_problem.grid.dx = problem.grid.dx * Ax;
homog_problem.grid.dy = problem.grid.dy * Ay;
homog_problem.grid.Lx = problem.grid.Lx;
homog_problem.grid.Ly = problem.grid.Ly;
homog_problem.Nx = NAx;
homog_problem.Ny = NAy;
homog_problem.nodeX = AnodeX;
homog_problem.nodeY = AnodeY;
homog_problem.stim_sites1 = Astim_sites1;
homog_problem.stim_sites2 = Astim_sites2;
homog_problem.cell_models = problem.cell_models;
homog_problem.model_assignments = Amodel_assignments;


end

