function FIGURE_SpiralActMaps
% This function loads in the activation data associated with the spiral
% wave problem, including both the finescale case and the homogenised cases
% present in the relevant .mat files

% Specify the ending timeframe (just used to ensure all time windows used
% are consistent
t_end = 4000;             % (in ms)
% Specify the window of time (relative to a starting moment of activation)
% that will be considered associated with the same activation event
t_window = 210;

% Specify homogenisation level to plot against finescale
homog_extent = 25;

% Load in the finescale activation data
load('spiral_finescale.mat', 'finescale_act_data', 'problem');

% Load in the homogenised model activation data
load('spiral_homog.mat', 'homog_act_data', 'homog_problems', 'homog_levels');

% Reference position for activation maps - time zero. This is scaled so 
% that zero is left/bottom, one is right/top
ref_X = 0.85;
ref_Y = 0.35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create re-labelled activation maps just to simplify coming code
% This reads out only homogenisation data for the extent requested above
AT_fine = finescale_act_data;
AT_homog = homog_act_data{ homog_levels == homog_extent };

% Also read out the homogenised problem
homog_problem = homog_problems{ homog_levels == homog_extent };

% Find number of nodes in this problem
N_nodes = length(problem.nodeX);

% Append an extra column of junk data onto the end of the activation time
% data (to record failures to find an activation event when using the "max
% trick"
AT_fine = [AT_fine, Inf * ones(N_nodes,1)];
AT_homog = [AT_homog, Inf * ones(N_nodes,1)];


% Find reference point
ref_fine = find( abs(problem.nodeX - problem.grid.Lx * ref_X) < 1e-6 & abs(problem.nodeY - problem.grid.Ly * ref_Y) < 1e-6, 1 );
ref_homog = find( abs(homog_problem.nodeX - homog_problem.grid.Lx * ref_X) < 1e-6 & abs(homog_problem.nodeY - homog_problem.grid.Ly * ref_Y) < 1e-6, 1 );

% Find number of activation events (including junk activation)
N_act_fine = size( AT_fine, 2 );
N_act_homog = size( AT_homog, 2 );

% Loop over each activation event at this site
looping = true;
c = 0;
while looping

    % Increment activation event counter
    c = c + 1;
    
    % Read out the starting time of this activation map
    t_ref = AT_fine(ref_fine, c);
    
    % Store this time in a vector
    fine_ref_vec(c) = t_ref;
    
    % If this is an actual activation event, and window remains within
    % maximum allowable time value, proceed
    if t_ref > 0 && t_ref < Inf && t_ref < t_end - t_window
        
        %%% FINESCALE
        
        % Find all subsequent activation events at each site
        % This uses a trick: max() will find first true value in each row
        [~,loc] = max( (AT_fine >= t_ref), [], 2);
        
        % Activation time is these time values minus reference time
        act_time = AT_fine( sub2ind( [N_nodes, N_act_fine], (1:N_nodes)', loc) ) - t_ref;
        
        % Non-activation events are set to a dummy value of -1
        % This will include those where next activation is after the
        % window, or there is no activation event (junk infinity values)
        act_time( act_time - t_window > 0 ) = -1;
        
        % Reshape into the dimensions of the problem for visualisation
        act_time = reshape( act_time, problem.Nx+1, problem.Ny+1)';
                
        % Visualise this activation map
        visualiseActivationMap(act_time, problem.occ_map);
        % Grab out axes (two per figure) 
        ax_objs = findobj(gcf,'Type','Axes');
        % Place title on topmost axes
        title(ax_objs(2),['t = ',num2str( round(t_ref, 2) ),' ms'], 'FontSize', 30);
        % Turn off all axis labels
        for k = 1:length(ax_objs)
            set(ax_objs(k),'xtick',[]);
            set(ax_objs(k),'xticklabels',[]);
            set(ax_objs(k),'ytick',[]);
            set(ax_objs(k),'yticklabels',[]);
        end
        caxis(ax_objs(2), [-t_window/256*1.01 t_window]);
        
        %%% HOMOGENISED
        
        % Calculate reference time using same point in the homogenised
        % model
        t_ref_homog = AT_homog(ref_homog, c);
        
        % Store this time in a vector
        homog_ref_vec(c) = t_ref_homog;
        
        % Find first activation events in the homogenised model that fall 
        % after the same reference time as the finescale 
        [~,loc] = max( (AT_homog >= t_ref_homog), [], 2 );
        
        % Activation time is these time values minus reference time
        act_time = AT_homog( sub2ind( [N_nodes, N_act_homog], (1:N_nodes)', loc) ) - t_ref_homog;
        
        % Remove any non-activation events, as above
        act_time( act_time - t_window > 0 ) = -1;
        
        % Reshape into the dimensions of the problem for visualisation
        act_time = reshape( act_time, homog_problem.Nx+1, homog_problem.Ny+1)';
        
        % Visualise this activation map
        visualiseActivationMap(act_time, homog_problem.occ_map);
        % Grab out axes (two per figure)
        ax_objs = findobj(gcf,'Type','Axes');
        % Place title on topmost axes
        title(ax_objs(2),['\Delta t = ',num2str( round(t_ref_homog - t_ref, 2) ),' ms'], 'FontSize', 30);
        % Turn off all axis labels
        for k = 1:length(ax_objs)
            set(ax_objs(k),'xtick',[]);
            set(ax_objs(k),'xticklabels',[]);
            set(ax_objs(k),'ytick',[]);
            set(ax_objs(k),'yticklabels',[]);
        end
        caxis(ax_objs(2), [-t_window/256*1.01 t_window]);
        
    % If this is not an activation event, terminate loop
    else
        looping = false;
    end
    
end

% Calculate frequencies - ignore last as it may not be complete
fine_freq = 1000 / mean( diff(fine_ref_vec(end-5:end-1) ) )
homog_freq = 1000 / mean( diff(homog_ref_vec(end-5:end-1) ) )

end