function createNozzleProblem(filename, channel_leftwidth, channel_rightwidth)
% This function creates a nozzle shape, a channel of a varying width at
% both ends. Emerging from the nozzle out into unobstructed tissue presents
% a situation of source/sink mismatch that can possibly trigger conduction
% block. 

% Define the diffusion tensor
D = [ 3, 0; 0, 3 ];    % Fibre-biased conduction

% Define the grid separation
mesh_separation = 10;

% The funnel shape is defined in terms of numbers of elements

% The number of elements in the whole domain
Nx = 400;
Ny = 100;

% The parameters controlling the funnel shape
margin = 100;             % Amount of margin on the left and right of the funnel structure

% Stimulus regions are the left and right edges of the domain (left
% primary, right secondary)
stim_width = 25;            % Width of stimulus regions (in elements)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Convert the grid separation to centimetres for consistent units
mesh_separation = mesh_separation / 10000;

% Find the domain size
Lx = Nx * mesh_separation;
Ly = Ny * mesh_separation;


% Determine x and y co-ordinates of all element centres
[X,Y] = meshgrid( linspace(0,Lx-mesh_separation,Nx) + mesh_separation/2,   linspace(0,Ly-mesh_separation,Ny) + mesh_separation/2 );

% Set up to have no channel on each side, just a wedge shape cut out of the
% block of obstruction
channel_xstart = margin+1;
channel_xend = Nx - margin;

% Create the funnel shape
occ_map = false(size(Y));
% For a wedge problem, only create a slope
slope = (channel_rightwidth/2 - channel_leftwidth/2) / ( channel_xend - channel_xstart );
occ_map(X(:)>mesh_separation * margin & X(:)<mesh_separation * (Nx - margin) & (Y(:) >= Ly/2 + channel_leftwidth/2 * mesh_separation + slope * (X(:) - mesh_separation*margin) - 100*eps | Y(:) <= Ly/2 - channel_leftwidth/2 * mesh_separation - slope * (X(:) - mesh_separation*margin) + 100*eps ) ) = true;


%%% Define volumes of each element
Vfrac = double(~occ_map);      % Volume fraction of one in all non-occupied elements (no elements are partially occupied)

%%% Create stimulus sites

% Nodes will be placed at element boundaries (vertex-centred finite volume)
% So first create node positions
[nodeX, nodeY] = meshgrid( linspace(0,Lx,Nx+1), linspace(0,Ly,Ny+1) );

% Initialise stimulus matrix to zeroes
stim_sites1 = false(size(nodeY));
stim_sites2 = false(size(nodeY));

% Set edges to be stimulus sites as requested
stim_sites1(:,1:stim_width) = true;
stim_sites2(:,end-stim_width+1:end) = true;



%%% Specify the cell model to use at all sites

% List cell models that will be used here
cell_models = {'TT3epi'};
% Assign models to cells (by number)
model_assignments = ones(size(nodeX));



%%% Process and save all data

% Read out base diffusivity levels from the diffusion tensor
D_xx = D(1,1);
D_xy = D(1,2);
D_yy = D(2,2);

% Create matrices of diffusion values, with zero in blocked regions
D_xx = D_xx * (~occ_map);
D_xy = D_xy * (~occ_map);
D_yy = D_yy * (~occ_map);

% Store problem details in the 'problem' structure
problem.occ_map = occ_map;
problem.D_tensor.D_xx = D_xx;
problem.D_tensor.D_xy = D_xy;
problem.D_tensor.D_yy = D_yy;
problem.Vfrac = Vfrac;
problem.grid.dx = mesh_separation;
problem.grid.dy = mesh_separation;
problem.grid.Lx = Lx;
problem.grid.Ly = Ly;
problem.Nx = Nx;
problem.Ny = Ny;
nodeX = nodeX'; nodeX = nodeX(:);
nodeY = nodeY'; nodeY = nodeY(:);
problem.nodeX = nodeX;
problem.nodeY = nodeY;
problem.stim_sites1 = stim_sites1;
problem.stim_sites2 = stim_sites2;
problem.cell_models = cell_models;
problem.model_assignments = model_assignments;

% Save the problem
save([filename,'.mat'],'problem');

end

