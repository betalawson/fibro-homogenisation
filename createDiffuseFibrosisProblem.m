function createDiffuseFibrosisProblem(filename, fibrosis_density, mesh_spacing)
% This function creates tissue of a single consistent property, then fills
% it with fibrotic occlusions of the requested density.

% Define the diffusion tensor
D = [ 3, 0; 0, 1 ];    % Fibre-biased conduction

% Define the physical size of the problem (in centimetres)
Lx = 5;
Ly = 0.5;

% Stimulus regions are the left and right edges of the domain (left
% primary, right secondary)
stim_width = 0.05;            % Width of stimulus regions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default value for the mesh spacing if not provided
if nargin < 3
    mesh_spacing = 100;
end

% Convert the mesh spacing to centimtres for consistent units
mesh_spacing = mesh_spacing / 10000;

% Use this new mesh spacing for both dx and dy
dx = mesh_spacing;
dy = mesh_spacing;

% Calculate the number of grid sites
Nx = ceil( Lx / dx );
Ny = ceil( Ly / dy );

% Initialise occupancy map
occ_map = false(Ny, Nx);

% Randomly turn off elements of the domain according to the requested
% fibrosis density
occ_map( rand(size(occ_map)) < fibrosis_density ) = true;


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
stim_sites1(nodeX <= stim_width) = true;



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
problem.grid.dx = dx;
problem.grid.dy = dy;
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

