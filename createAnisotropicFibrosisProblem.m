function createAnisotropicFibrosisProblem(filename, fibrosis_density, mesh_spacing, orientation)
% This function creates tissue of a single consistent property, then fills
% it with fibrotic occlusions of the requested density. The fibrotic
% occlusions are arranged as barriers oriented either vertically or
% horizontally (according to 'x' or 'y' value for supplied input
% orientation). This creates fibrotic patterns with an evident anisotropy

% Define the diffusion tensor
D = [ 3, 0; 0, 1 ];

% Define the physical size of the problem (in centimetres)
Lx = 5;
Ly = 0.5;

% Stimulus regions are the left and right edges of the domain (left
% primary, right secondary)
stim_width = 0.05;            % Width of stimulus regions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert the mesh spacing to centimtres for consistent units
mesh_spacing = mesh_spacing / 10000;

% Use this new mesh spacing for both dx and dy
dx = mesh_spacing;
dy = mesh_spacing;

% Calculate the number of grid sites
Nx = ceil( Lx / dx );
Ny = ceil( Ly / dy );

%%% Create a map of where fibrosis is, with density approximately equal to
%%% "fibrosis_density"
occ_map = false(Ny,Nx);
density = 0;
while density < fibrosis_density

   % Select a random location
   x_start = ceil(rand * Nx);
   y_start = ceil(rand * Ny);
   
   % Define length of this fibre
   neg_length = 4;
   pos_length = 4;
   
   % Place this fibre
   switch orientation
       case 'x'
           occ_map(y_start, max([1, x_start - neg_length]):min([Nx, x_start + pos_length])) = true;
       case 'y'
           occ_map(max([1, y_start - neg_length]):min([Ny, y_start + pos_length]),x_start) = true;
   end
   
   % Update density
   density = sum(occ_map(:)) / numel(occ_map);
    
end


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
stim_sites2( sqrt(nodeX.^2 + nodeY.^2) <= stim_width ) = true;



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
problem.nodeX = nodeX;
problem.nodeY = nodeY;
problem.stim_sites1 = stim_sites1;
problem.stim_sites2 = stim_sites2;
problem.cell_models = cell_models;
problem.model_assignments = model_assignments;

% Save the problem
save([filename,'.mat'],'problem');

end

