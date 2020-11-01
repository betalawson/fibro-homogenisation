function [Jx_bar, Jy_bar] = solveClosureProblemPeriodic(problem, Edir, skin_map)
% This function takes an input problem (specified in the traditional
% 'problem' format used by the monodomain solver code) and solves the
% steady-state diffusion problem with no reaction term, but subject to an
% external electric field, for which Edir specifies the direction (x or y).
% Linear boundary conditions are applied (no-flux on internal boundaries).
% If a skin map is not provided, or provided as an empty vector, all sites
% will be used in calculating averaged fluxes. Otherwise, only those marked
% as non-skin sites will be used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read out the required problem information
occ_map = problem.occ_map;
D_xx = problem.D_tensor.D_xx;
D_xy = problem.D_tensor.D_xy;
D_yy = problem.D_tensor.D_yy;
Vfrac = problem.Vfrac;
dx = problem.grid.dx;
dy = problem.grid.dy;

% Read out the size of the input matrices
[Ny, Nx] = size(D_xx);


% Make the border surrounding the problem domain a copy of the opposite
% edges to obtain a periodic substructure. This simply ensures that flux is
% calculated correctly in border regions without the need for awkward loops
% over boundaries and corners
occ_map_ext = [ [occ_map(Ny,Nx), occ_map(Ny,1:Nx), occ_map(Ny,1) ];   % Top row
                [occ_map(1:Ny,Nx), occ_map, occ_map(1:Ny,1)];         % Main rows
                [occ_map(1,Nx), occ_map(1,1:Nx), occ_map(1,1) ] ];    % Bottom row


% Perform the same process for the volume fraction
Vfrac_ext = [ [Vfrac(Ny,Nx), Vfrac(Ny,1:Nx), Vfrac(Ny,1) ];   % Top row
                [Vfrac(1:Ny,Nx), Vfrac, Vfrac(1:Ny,1)];         % Main rows
                [Vfrac(1,Nx), Vfrac(1,1:Nx), Vfrac(1,1) ] ];    % Bottom row


% Match this with the periodic extensions of the diffusion tensors
D_xx_ext = [ [D_xx(Ny,Nx), D_xx(Ny,1:Nx), D_xx(Ny,1) ];   % Top row
             [D_xx(1:Ny,Nx), D_xx, D_xx(1:Ny,1)];         % Main rows
             [D_xx(1,Nx), D_xx(1,1:Nx), D_xx(1,1) ] ];    % Bottom row
D_xy_ext = [ [D_xy(Ny,Nx), D_xy(Ny,1:Nx), D_xy(Ny,1) ];   % Top row
             [D_xy(1:Ny,Nx), D_xy, D_xy(1:Ny,1)];         % Main rows
             [D_xy(1,Nx), D_xy(1,1:Nx), D_xy(1,1) ] ];    % Bottom row
D_yy_ext = [ [D_yy(Ny,Nx), D_yy(Ny,1:Nx), D_yy(Ny,1) ];   % Top row
             [D_yy(1:Ny,Nx), D_yy, D_yy(1:Ny,1)];         % Main rows
             [D_yy(1,Nx), D_yy(1,1:Nx), D_yy(1,1) ] ];    % Bottom row

         
% A very similar calculation calculates the volume of each control volume
CV_vols = dx * dy * ( ~occ_map_ext(1:Ny+1, 1:Nx+1) + ~occ_map_ext(2:Ny+2, 1:Nx+1) + ~occ_map_ext(1:Ny+1, 2:Nx+2) + ~occ_map_ext(2:Ny+2, 2:Nx+2) ) / 4;

% Convert these matrices into a vector for ease of use in the matrix system
Vfrac_v = Vfrac'; Vfrac_v = Vfrac_v(:);
Vfrac_ext = Vfrac_ext'; Vfrac_ext = Vfrac_ext(:);
CV_vols = CV_vols'; CV_vols = CV_vols(:);

% Define the number of nodes in the normal and extended node grids
Nn = (Nx+1)*(Ny+1);

% Calculate the dependences of the bilinearly interpolated flux values at
% the four control volume boundary midpoints for each non-fibrotic element.
% By dependencies, this refers to how these fluxes depend on each of the
% four surrounding note points. To do this, we use a vectorised approach
% (i.e. occupancy matrix becomes a vector)

% Convert occupancy and diffusivities to vectors
occ_v = occ_map'; occ_v = occ_v(:);
occ_v_ext = occ_map_ext'; occ_v_ext = occ_v_ext(:);
D_xx_v = D_xx'; D_xx_v = D_xx_v(:);
D_xy_v = D_xy'; D_xy_v = D_xy_v(:);
D_yy_v = D_yy'; D_yy_v = D_yy_v(:);
D_xx_extv = D_xx_ext'; D_xx_extv = D_xx_extv(:);
D_xy_extv = D_xy_ext'; D_xy_extv = D_xy_extv(:);
D_yy_extv = D_yy_ext'; D_yy_extv = D_yy_extv(:);

% Vertical fluxes through the control volume boundaries [dl, dr, ul, ur]
J_W = Vfrac_ext .* ~occ_v_ext .* ( 1 / (dx*dy) ) .* ( [ -D_xy_extv * dy/2 - D_yy_extv * 3*dx/4, D_xy_extv * dy/2 - D_yy_extv * dx/4, -D_xy_extv * dy/2 + D_yy_extv * 3*dx/4, D_xy_extv * dy/2 + D_yy_extv * dx/4 ] );
J_E = Vfrac_ext .* ~occ_v_ext .* ( 1 / (dx*dy) ) .* ( [ -D_xy_extv * dy/2 - D_yy_extv * dx/4, D_xy_extv * dy/2 - D_yy_extv * 3*dx/4, -D_xy_extv * dy/2 + D_yy_extv * dx/4, D_xy_extv * dy/2 + D_yy_extv * 3*dx/4 ] );
J_S = Vfrac_ext .* ~occ_v_ext .* ( 1 / (dx*dy) ) .* ( [ -D_xy_extv * dx/2 - D_xx_extv * 3*dy/4, -D_xy_extv * dx/2 + D_xx_extv * 3*dy/4, D_xy_extv * dx/2 - D_xx_extv * dy/4, D_xy_extv * dx/2 + D_xx_extv * dy/4 ] );
J_N = Vfrac_ext .* ~occ_v_ext .* ( 1 / (dx*dy) ) .* ( [ -D_xy_extv * dx/2 - D_xx_extv * dy/4, -D_xy_extv * dx/2 + D_xx_extv * dy/4, D_xy_extv * dx/2 - D_xx_extv * 3*dy/4, D_xy_extv * dx/2 + D_xx_extv * 3*dy/4 ] );

% First create lists of the i and j co-ordinates of all *node* points
j = ceil( (1:Nn) / (Nx+1) );
i = ( (1:Nn) - (j-1)*(Nx+1) );

% Create a vector of the corresponding position in the *element* grid
% (which is Ny+2 by Nx+2 after extension with boundaries). This position
% refers to the element down and left from the node point
eloc_dl_i = i;
eloc_dl_j = j;
eloc_dr_i = i+1;
eloc_dr_j = j;
eloc_ul_i = i;
eloc_ul_j = j+1;
eloc_ur_i = i+1;
eloc_ur_j = j+1;

% Now set all values outside of the 'real' domain to be the corresponding
% elements on the opposite side instead
eloc_dl_i( eloc_dl_i == 1) = Nx+1;
eloc_dr_i( eloc_dr_i == 1) = Nx+1;
eloc_ul_i( eloc_ul_i == 1) = Nx+1;
eloc_ur_i( eloc_ur_i == 1) = Nx+1;
eloc_dl_i( eloc_dl_i == Nx+2) = 2;
eloc_dr_i( eloc_dr_i == Nx+2) = 2;
eloc_ul_i( eloc_ul_i == Nx+2) = 2;
eloc_ur_i( eloc_ur_i == Nx+2) = 2;
eloc_dl_j( eloc_dl_j == 1) = Ny+1;
eloc_dr_j( eloc_dr_j == 1) = Ny+1;
eloc_ul_j( eloc_ul_j == 1) = Ny+1;
eloc_ur_j( eloc_ur_j == 1) = Ny+1;
eloc_dl_j( eloc_dl_j == Ny+2) = 2;
eloc_dr_j( eloc_dr_j == Ny+2) = 2;
eloc_ul_j( eloc_ul_j == Ny+2) = 2;
eloc_ur_j( eloc_ur_j == Ny+2) = 2;

% Now convert these to element-wise co-ordinates
eloc_dl = ( (eloc_dl_j - 1 ) * (Nx+2) + eloc_dl_i )';
eloc_dr = ( (eloc_dr_j - 1 ) * (Nx+2) + eloc_dr_i )';
eloc_ul = ( (eloc_ul_j - 1 ) * (Nx+2) + eloc_ul_i )';
eloc_ur = ( (eloc_ur_j - 1 ) * (Nx+2) + eloc_ur_i )';

% Use i and j to create a list of co-ordinates for each node point, and its
% neighbours
P_i = i; P_j = j;
E_i = i+1; E_j = j;
W_i = i-1; W_j = j;
N_i = i; N_j = j+1;
S_i = i; S_j = j-1;
NE_i = i+1; NE_j = j+1;
SE_i = i+1; SE_j = j-1;
NW_i = i-1; NW_j = j+1;
SW_i = i-1; SW_j = j-1;

% Apply periodic looping, so that all node references outside grid loop
% around to their counterpart
E_i(E_i == Nx+2) = 2;
W_i(W_i == 0) = Nx;
N_j(N_j == Ny+2) = 2;
S_j(S_j == 0) = Ny;
NE_i(NE_i == Nx+2) = 2;
SE_i(SE_i == Nx+2) = 2;
NW_i(NW_i == 0) = Nx;
SW_i(SW_i == 0) = Nx;
NE_j(NE_j == Ny+2) = 2;
NW_j(NW_j == Ny+2) = 2;
SE_j(SE_j == 0) = Ny;
SW_j(SW_j == 0) = Ny;

% Convert these into element-wise co-ordinates
P = ( P_j - 1 ) * (Nx+1) + P_i;
E = ( E_j - 1 ) * (Nx+1) + E_i;
W = ( W_j - 1 ) * (Nx+1) + W_i;
N = ( N_j - 1 ) * (Nx+1) + N_i;
S = ( S_j - 1 ) * (Nx+1) + S_i;
NE = ( NE_j - 1 ) * (Nx+1) + NE_i;
SE = ( SE_j - 1 ) * (Nx+1) + SE_i;
NW = ( NW_j - 1 ) * (Nx+1) + NW_i;
SW = ( SW_j - 1 ) * (Nx+1) + SW_i;

% Create lists of node locations along the boundaries
Esites = Nx+1 : Nx+1 : (Nx+1)*(Ny+1);
Wsites = 1 : Nx+1 : Ny*(Nx+1)+1;
Nsites = Ny*(Nx+1) + (1 : Nx+1);
Ssites = 1 : Nx+1;


% Now, the control volume formulation of the diffusive part of the problem
% is used to build the matrix, in terms of the above fluxes. The sparse
% matrix is created by first creating vectors [i,j,v] that define the row,
% column and value of the elements, respectively. This saves a great deal
% of memory and time
i_v = []; j_v = []; val_v = [];

% Node's dependence on itself
i_v = [i_v, P]; j_v = [j_v, P];
val_v = [val_v; ( dx/2 * J_W(eloc_ur,1) + dy/2 * J_S(eloc_ur,1) + dx/2 * J_E(eloc_ul,2) - dy/2 * J_S(eloc_ul,2) - dx/2 * J_W(eloc_dr,3) + dy/2 * J_N(eloc_dr,3) - dx/2 * J_E(eloc_dl,4) - dy/2 * J_N(eloc_dl,4) ) ./ CV_vols ];

% Node's dependence on North-East node
i_v = [i_v, P]; j_v = [j_v, NE];
val_v = [val_v; ( dx/2 * J_W(eloc_ur,4) + dy/2 * J_S(eloc_ur,4) ) ./ CV_vols];

% Node's dependence on North-West node
i_v = [i_v, P]; j_v = [j_v, NW];
val_v = [val_v; ( dx/2 * J_E(eloc_ul,3) - dy/2 * J_S(eloc_ul,3) ) ./ CV_vols];

% Node's dependence on South-East node
i_v = [i_v, P]; j_v = [j_v, SE];
val_v = [val_v; (- dx/2 * J_W(eloc_dr,2) + dy/2 * J_N(eloc_dr,2) ) ./ CV_vols];

% Node's dependence on South-West node
i_v = [i_v, P]; j_v = [j_v, SW];
val_v = [val_v; (- dx/2 * J_E(eloc_dl,1) - dy/2 * J_N(eloc_dl,1) ) ./ CV_vols];

% Node's dependence on East node
i_v = [i_v, P]; j_v = [j_v, E];
val_v = [val_v; (dx/2 * J_W(eloc_ur,2) + dy/2 * J_S(eloc_ur,2) - dx/2 * J_W(eloc_dr,4) + dy/2 * J_N(eloc_dr,4) ) ./ CV_vols];

% Node's dependence on West node
i_v = [i_v, P]; j_v = [j_v, W];
val_v = [val_v; (dx/2 * J_E(eloc_ul,1) - dy/2 * J_S(eloc_ul,1) - dx/2 * J_E(eloc_dl,3) - dy/2 * J_N(eloc_dl,3) ) ./ CV_vols];

% Node's dependence on North node
i_v = [i_v, P]; j_v = [j_v, N];
val_v = [val_v; (dx/2 * J_W(eloc_ur,3) + dy/2 * J_S(eloc_ur,3) + dx/2 * J_E(eloc_ul,4) - dy/2 * J_S(eloc_ul,4) ) ./ CV_vols];

% Node's dependence on South node
i_v = [i_v, P]; j_v = [j_v, S];
val_v = [val_v; (- dx/2 * J_W(eloc_dr,1) + dy/2 * J_N(eloc_dr,1) - dx/2 * J_E(eloc_dl,2) - dy/2 * J_N(eloc_dl,2) ) ./ CV_vols];

% Before creating the matrix, set any values that appear to be zero (but
% are not due to rounding error) to explicitly zero so they are not stored
val_v( abs(val_v) < 100*eps ) = 0;

% Now create the matrix
A = sparse(i_v, j_v, val_v, Nn, Nn);

% Initialise with no RHS vector
b = zeros(Nn,1);

% Now, add a component of flux that is contributed by every unoccupied
% neighbour. If all neighbours are unoccupied, this will result in a zero
% value for constant diffusion coefficients
switch Edir
    case 'x'
        b = -( (-D_xx_extv(eloc_dl) + D_xx_extv(eloc_dr) - D_xx_extv(eloc_ul) + D_xx_extv(eloc_ur) ) * (dy / 2)  +  (-D_xy_extv(eloc_dl) - D_xy_extv(eloc_dr) + D_xy_extv(eloc_ul) + D_xy_extv(eloc_ur) ) * (dx / 2) ) ./ CV_vols;
    case 'y'
        b = -( (-D_xy_extv(eloc_dl) + D_xy_extv(eloc_dr) - D_xy_extv(eloc_ul) + D_xy_extv(eloc_ur) ) * (dy / 2)  +  (-D_yy_extv(eloc_dl) - D_yy_extv(eloc_dr) + D_yy_extv(eloc_ul) + D_yy_extv(eloc_ur) ) * (dx / 2) ) ./ CV_vols;
end

% Apply the periodic boundary condition by enforcing that all derivative
% nodes must equal their periodic equivalent:

% First, zero out all corresponding rows (also RHS vector values)
A( [Esites, Nsites], : ) = 0;
b( [Esites, Nsites] ) = 0;

% Now fill in all rows with a pattern of <this site> + -1<match site> = 0
A(sub2ind([(Ny+1)*(Nx+1), (Ny+1)*(Nx+1)], Esites, Esites) ) = 1;
A(sub2ind([(Ny+1)*(Nx+1), (Ny+1)*(Nx+1)], Esites, Wsites) ) = -1;
A(sub2ind([(Ny+1)*(Nx+1), (Ny+1)*(Nx+1)], Nsites, Nsites) ) = 1;
A(sub2ind([(Ny+1)*(Nx+1), (Ny+1)*(Nx+1)], Nsites, Ssites) ) = -1;
% Make sure all corners also match (somewhat included above, but to be
% safe)
A( (Ny+1)*(Nx+1), (Ny+1)*(Nx+1) ) = 1;
A( (Ny+1)*(Nx+1), 1 ) = -1;
A( Ny*(Nx+1), Ny*(Nx+1) ) = 1;
A( Ny*(Nx+1), 1 ) = 1;
A( Nx+1, Nx+1 ) = 1;
A( Nx+1, 1 ) = -1;
        
% Track which nodes are 'active' (not buried completely in fibrosis)
active = ~all([occ_v_ext(eloc_dl), occ_v_ext(eloc_dr), occ_v_ext(eloc_ul), occ_v_ext(eloc_ur) ], 2);

% Reduce diffusive update matrix and vector according to these nodes. This
% also removes rows of zeros from the matrix, allowing use of backslash
A = A(active, active);
b = b(active);


% Now, solve the system. lsqminnorm is used, because isolated regions that
% are disconnected from boundaries are 'ill-posed' (only defined within a
% constant)
w_active = lsqminnorm(A,b);

% Create a vectorised list of all closure variable values, using NaN for 
% nodes that are inactive
w = nan( (Ny+1)*(Nx+1),1);
w(active) = w_active;



% In preparation for the average flux calculations, convert the skin_map
% input to a vector, if it was provided. Otherwise, just mark all sites as
% not skin so the averaging can be performed using all sites
if nargin > 2 && ~isempty(skin_map)
    skin_map = skin_map'; skin_map = skin_map(:);
else
    skin_map = false(size(occ_v));
end


% Grab out a list of numerical indices for unoccupied elements that are not
% skin sites - these are the sites that will be used for averaging
avg_ele = find(~occ_v);
avg_ele = avg_ele(~skin_map(avg_ele));

% Create a list of the nodes corresponding to each of these elements
dl = avg_ele + floor( (avg_ele-1) / Nx );
dr = dl + 1;
ul = dl + Nx+1;
ur = dl + Nx+1 + 1;

% Now find the average flux in the x and y directions (bilinear
% interpolation)
Jx_bar = sum(  Vfrac_v(avg_ele) .* (  dy / 2 * D_xx_v(avg_ele) .* ( w(dr) + w(ur) - w(dl) - w(ul) )  + dx / 2 * D_xy_v(avg_ele) .* ( w(ul) + w(ur) - w(dl) - w(dr) ) ) ) / sum( Vfrac_v(avg_ele) * dx * dy );
Jy_bar = sum(  Vfrac_v(avg_ele) .* (  dy / 2 * D_xy_v(avg_ele) .* ( w(dr) + w(ur) - w(dl) - w(ul) )  + dx / 2 * D_yy_v(avg_ele) .* ( w(ul) + w(ur) - w(dl) - w(dr) ) ) ) / sum( Vfrac_v(avg_ele) * dx * dy );

% Add on the flux in the direction of the imposed unit electric field
switch Edir
    case 'x'
        Jx_bar = Jx_bar + sum(  Vfrac_v(avg_ele) .* D_xx_v(avg_ele) ) / sum( Vfrac_v(avg_ele) );
        Jy_bar = Jy_bar + sum(  Vfrac_v(avg_ele) .* D_xy_v(avg_ele) ) / sum( Vfrac_v(avg_ele) );
    case 'y'
        Jx_bar = Jx_bar + sum(  Vfrac_v(avg_ele) .* D_xy_v(avg_ele) ) / sum( Vfrac_v(avg_ele) );
        Jy_bar = Jy_bar + sum(  Vfrac_v(avg_ele) .* D_yy_v(avg_ele) ) / sum( Vfrac_v(avg_ele) );
end

% Code is done, return values Jx_bar, Jy_bar

end