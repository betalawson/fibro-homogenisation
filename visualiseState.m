function visualiseState( V, Nx, Ny, occ_map, t )
% This function visualises the membrane potential map

% Internal flag specifies whether to visualise as 'surfaces' or points
vis_surfaces = false;

% Convert vector V back into a matrix (use transpose for correct
% orientation because reshape works 'columnwise')
V = reshape(V, Nx+1, Ny+1)';

% Load colormap
load('extra_colormaps.mat','plasma');
Vclr = plasma;

if vis_surfaces

    % Append grey for visualisation of regions without nodes
    Vclr = [ [0.2, 0.2, 0.2]; Vclr]; 
    
    % Calculate *element* Voltage values as the average of node corners
    % Some will be NaNs, but we are happy for these to propagate because they
    % will only 'infect' sites that are occupied anyway
    Ve = ( V(1:Ny,1:Nx) + V(2:Ny+1,1:Nx) + V(1:Ny,2:Nx+1) + V(2:Ny+1,2:Nx+1) ) / 4;

    % Make sure all occupied sites won't show up in the plot
    Ve(occ_map) = NaN;

    % Visualise the remainder
    imagesc(flipud(Ve));
    
    % Visualisation properties
    colormap(Vclr);
    whitebg([0.2 0.2 0.2]);
    caxis([-90 40]);
    axis equal;

    % Add title
    title(['Membrane potential map after t = ',num2str(t),' ms'],'FontSize', 24);
    
else
    
    % Set dotsize according to problem dimension
    dotsize = max([5, 500/Nx, 500/Ny]);
    
    % Visualise nodes using scatter
    [X,Y] = meshgrid( 1:Nx+1, 1:Ny+1 );
    scatter( X(:), Y(:), dotsize, V(:), 'filled' );
    xlim([0, Nx+2]);
    ylim([0, Ny+2]);
    
    % Visualisation properties
    colormap(Vclr);
    whitebg([1 1 1]);
    caxis([-90 40]);
    axis equal;

    % Add title
    title(['Membrane potential map after t = ',num2str(t),' ms'],'FontSize', 24);

end

