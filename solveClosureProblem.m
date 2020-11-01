function [Jx_bar, Jy_bar] = solveClosureProblem(problem, Edir, boundary_type, skin_map)
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


% Make the border surrounding the problem domain fibrotic. This simply
% ensures that flux is calculated correctly in border regions without the
% need for awkward loops over boundaries and corners
occ_map_ext = [true(1,Nx+2); [true(Ny,1), occ_map, true(Ny,1)]; true(1,Nx+2) ];

% Perform the same process for the volume fraction
Vfrac_ext = [zeros(1,Nx+2); [zeros(Ny,1), Vfrac, zeros(Ny,1)]; zeros(1,Nx+2) ];

% Match this with a layer of zero diffusivity around the diffusivity maps
% provided
D_xx_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_xx, zeros(Ny,1)]; zeros(1,Nx+2) ];
D_xy_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_xy, zeros(Ny,1)]; zeros(1,Nx+2) ];
D_yy_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_yy, zeros(Ny,1)]; zeros(1,Nx+2) ];

% A very similar calculation calculates the volume of each control volume
CV_vols = dx * dy * ( ~occ_map_ext(1:Ny+1, 1:Nx+1) + ~occ_map_ext(2:Ny+2, 1:Nx+1) + ~occ_map_ext(1:Ny+1, 2:Nx+2) + ~occ_map_ext(2:Ny+2, 2:Nx+2) ) / 4;

% Convert these matrices into a vector for ease of use in the matrix system
Vfrac_v = Vfrac'; Vfrac_v = Vfrac_v(:);
Vfrac_ext = Vfrac_ext'; Vfrac_ext = Vfrac_ext(:);
CV_vols = CV_vols'; CV_vols = CV_vols(:);

% Define the number of nodes in the normal and extended node grids
Nn = (Nx+1)*(Ny+1);
Nf = (Nx+3)*(Ny+3);

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
eloc_dl = ( (j-1) * (Nx+2) + i )';

% Calculate the corresponding positions of other elements accordingly
eloc_dr = eloc_dl + 1;
eloc_ul = eloc_dl + (Nx+2);
eloc_ur = eloc_dl + (Nx+2) + 1;

% Also create a remapped vector that converts the list of all node
% locations to their locations in a larger (Nx+3) by (Ny+3) grid of nodes
% (so that everything can be done with matrix indexing)
nlist = ( j * (Nx+3) + i + 1 );

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
i_v = [i_v, nlist]; j_v = [j_v, nlist];
val_v = [val_v; ( dx/2 * J_W(eloc_ur,1) + dy/2 * J_S(eloc_ur,1) + dx/2 * J_E(eloc_ul,2) - dy/2 * J_S(eloc_ul,2) - dx/2 * J_W(eloc_dr,3) + dy/2 * J_N(eloc_dr,3) - dx/2 * J_E(eloc_dl,4) - dy/2 * J_N(eloc_dl,4) ) ./ CV_vols];

% Node's dependence on North-East node
i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)+1];
val_v = [val_v; ( dx/2 * J_W(eloc_ur,4) + dy/2 * J_S(eloc_ur,4) ) ./ CV_vols];

% Node's dependence on North-West node
i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)-1];
val_v = [val_v; ( dx/2 * J_E(eloc_ul,3) - dy/2 * J_S(eloc_ul,3) ) ./ CV_vols];

% Node's dependence on South-East node
i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)+1];
val_v = [val_v; (- dx/2 * J_W(eloc_dr,2) + dy/2 * J_N(eloc_dr,2) ) ./ CV_vols];

% Node's dependence on South-West node
i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)-1];
val_v = [val_v; (- dx/2 * J_E(eloc_dl,1) - dy/2 * J_N(eloc_dl,1) ) ./ CV_vols];

% Node's dependence on East node
i_v = [i_v, nlist]; j_v = [j_v, nlist+1];
val_v = [val_v; (dx/2 * J_W(eloc_ur,2) + dy/2 * J_S(eloc_ur,2) - dx/2 * J_W(eloc_dr,4) + dy/2 * J_N(eloc_dr,4) ) ./ CV_vols];

% Node's dependence on West node
i_v = [i_v, nlist]; j_v = [j_v, nlist-1];
val_v = [val_v; (dx/2 * J_E(eloc_ul,1) - dy/2 * J_S(eloc_ul,1) - dx/2 * J_E(eloc_dl,3) - dy/2 * J_N(eloc_dl,3) ) ./ CV_vols];

% Node's dependence on North node
i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)];
val_v = [val_v; (dx/2 * J_W(eloc_ur,3) + dy/2 * J_S(eloc_ur,3) + dx/2 * J_E(eloc_ul,4) - dy/2 * J_S(eloc_ul,4) ) ./ CV_vols];

% Node's dependence on South node
i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)];
val_v = [val_v; (- dx/2 * J_W(eloc_dr,1) + dy/2 * J_N(eloc_dr,1) - dx/2 * J_E(eloc_dl,2) - dy/2 * J_N(eloc_dl,2) ) ./ CV_vols];

% Before creating the matrix, set any values that appear to be zero (but
% are not due to rounding error) to explicitly zero so they are not stored
val_v( abs(val_v) < 100*eps ) = 0;

% Now create the matrix
A = sparse(i_v, j_v, val_v, Nf, Nf);

% Now, grab out only the rows of the matrix that correspond to real nodes
% (i.e. without the dummy boundaries)
A = A(nlist,nlist);

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

% Track which nodes are 'active' (not buried completely in fibrosis)
active = ~all([occ_v_ext(eloc_dl), occ_v_ext(eloc_dr), occ_v_ext(eloc_ul), occ_v_ext(eloc_ur) ], 2);

% Apply the linear boundary condition by specifying Drichlet boundary 
% condition of zero on all boundary sites
        
% Read out the active nodes on all boundaries
AEsites = Esites(active(Esites));
AWsites = Wsites(active(Wsites));
ANsites = Nsites(active(Nsites));
ASsites = Ssites(active(Ssites));
        
switch boundary_type
    
    case 'linear'

        % Apply Drichlet condition to all of the boundaries
        bound_sites = [AEsites, AWsites, ANsites, ASsites];
        A(bound_sites,:) = 0;
        A(sub2ind([(Ny+1)*(Nx+1), (Ny+1)*(Nx+1)], bound_sites, bound_sites) ) = 1;
        b(bound_sites) = 0;
        
    case 'confined'
   
        % Apply Drichlet condition only to the boundaries that create the
        % imposed gradient (depends on direction of imposed gradient)
        switch Edir
            
            case 'x'
                bound_sites = [AEsites, AWsites];
                A(bound_sites,:) = 0;
                A(sub2ind([(Ny+1)*(Nx+1), (Ny+1)*(Nx+1)], bound_sites, bound_sites) ) = 1;
                b(bound_sites) = 0;
                
            case 'y'
                bound_sites = [ANsites, ASsites];
                A(bound_sites,:) = 0;
                A(sub2ind([(Ny+1)*(Nx+1), (Ny+1)*(Nx+1)], bound_sites, bound_sites) ) = 1;
                b(bound_sites) = 0;
                
        end
        
end


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

if nargin <= 3 || isempty(skin_map)
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