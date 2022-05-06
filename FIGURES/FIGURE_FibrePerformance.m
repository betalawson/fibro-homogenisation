function FIGURE_FibrePerformance
% This function plots the figure showing the performance of homogenisation
% for the different extents of homogenisation and different densities of
% fibrosis.

% Velocity measuring information
Lx = 5;                      % Length of fibres (cm)
start_point = 0.3;          % Proportion along fibre to use as start point
end_point = 0.8;            % Proportion along fibre to use as end point

% List the homogenisation types used
homog_types = {'linear', 'confined', 'periodic'};

% List the fibrosis types used
fibrosis_types = {'diffuse', 'parallel', 'perp'};

% Define a list of colours
colors = [ [0.1 0.1 0.7];
           [0.45 0.45 0.95];
           [0.95 0.45 0.45];
           [0.8 0.15 0.15]
         ];


% Use the "original grid" simulations to calculate homogenised velocities
% in a manner free of grid considerations
for j = 1:length(fibrosis_types)
    for i = 1:length(homog_types)
        
        % Load in the base data for this type of fibrosis
        load([fibrosis_types{j},'_finescale.mat'], 'fine_act_maps', 'fibrosis_densities');
        
        % Calculate the columns that correspond to the 35% and 85% positions along
        % the fibre - Nx and Ny are numbers of nodes, not elements. This is
        % only inside the loop because it uses the read-in data to
        % determine mesh size
        [Ny,Nx] = size(fine_act_maps{1});
        start_loc = round((Nx-1)*start_point)+1;
        end_loc = round((Nx-1)*end_point)+1;
        
        % Load in the data for this type of homogenisation
        %load([fibrosis_types{j},'_OG_',homog_types{i}], 'homog_act_maps');
        %OG_homog_maps = homog_act_maps;
        load([fibrosis_types{j},'_RG_',homog_types{i}],'grid_sizes', 'homog_act_maps');
        RG_homog_maps = homog_act_maps;
        
        % Convert the grid sizes into the numbers of elements
        homog_factors = (grid_sizes / 10);
        
        % Loop over the different extents of homogenisation applied
        for k = 1:length(grid_sizes)
            
            % Calculate the correct read-out points for the reduced
            % grid homogenised results (which have different numbers of
            % node points)
            RG_start_loc = (start_loc - 1) / homog_factors(k) + 1;
            RG_end_loc = (end_loc - 1) / homog_factors(k) + 1;
            
            % Loop over the different densities of fibrosis
            for m = 1:size(RG_homog_maps,2)
                
                % Make sure non-activated sites are recorded as NaN, not -1
                fine_act_maps{m}( fine_act_maps{m} == -1 ) = NaN;
                RG_homog_maps{k,m}( RG_homog_maps{k,m} == -1 ) = NaN;
                %OG_homog_maps{k,m}( OG_homog_maps{k,m} == -1 ) = NaN;
                
                % Calculate the velocity
                base_vels(m) = Lx * (end_point - start_point) / (nanmean(fine_act_maps{m}(:,end_loc)) - nanmean(fine_act_maps{m}(:,start_loc))) * 1000;      % Converted into cm/s
                %OG_homog_vels(k,m) = Lx * (end_point - start_point) / (nanmean(OG_homog_maps{k,m}(:,end_loc)) - nanmean(OG_homog_maps{k,m}(:,start_loc))) * 1000;      % Converted into cm/s
                RG_homog_vels(k,m) = Lx * (end_point - start_point) / (nanmean(RG_homog_maps{k,m}(:,RG_end_loc)) - nanmean(RG_homog_maps{k,m}(:,RG_start_loc))) * 1000;      % Converted into cm/s
                
            end
        end
        
        %%% Create the figure - original grid
%         figure('units','normalized','OuterPosition',[0 0 1 1]);
%         hold on;
%         
%         % Plot the data - also create legend while looping
%         leg_text = {'Fine scale'};
%         plot(fibrosis_densities, base_vels, 'k.', 'MarkerSize', 36);
%         for k = 1:length(grid_sizes)
%             plot(fibrosis_densities, OG_homog_vels(k,:), '.', 'MarkerSize', 30, 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', colors(k,:) );
%             leg_text{k+1} = [num2str(homog_factors(k)),'x',num2str(homog_factors(k)),' Homogenisation'];
%         end
%         
%         % Add labels
%         xlabel('Proportion of Tissue Obstructed', 'FontSize', 22);
%         ylabel('Conduction Velocity (cm/s)', 'FontSize', 22);
%         
%         % Set up limits
%         xlim([fibrosis_densities(1) - (fibrosis_densities(2) - fibrosis_densities(1)), fibrosis_densities(end) + (fibrosis_densities(2) - fibrosis_densities(1))]);
%         
%         % Set up axes
%         set(gca,'FontSize', 20);
%         
%         % Append legend
%         legend(leg_text);
        
        
        %%% Create the figure - reduced grid
        figure('units','normalized','OuterPosition',[0 0 1 1]);
        hold on;
        
        % Plot the data - also create legend while looping
        leg_text = {'Fine scale'};
        plot(fibrosis_densities, base_vels, 'k.', 'MarkerSize', 36);
        for k = 1:length(grid_sizes)
            plot(fibrosis_densities, RG_homog_vels(k,:), '.', 'MarkerSize', 30, 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', colors(k,:) );
            leg_text{k+1} = [num2str(homog_factors(k)),'x',num2str(homog_factors(k)),' Homogenisation'];
        end
        
        % Add labels
        xlabel('Proportion of Tissue Obstructed', 'FontSize', 22);
        ylabel('Conduction Velocity (cm/s)', 'FontSize', 22);
        
        % Set up limits
        xlim([fibrosis_densities(1) - (fibrosis_densities(2) - fibrosis_densities(1)), fibrosis_densities(end) + (fibrosis_densities(2) - fibrosis_densities(1))]);
        
        % Set up axes
        set(gca,'FontSize', 20);
        
        % Append legend
        legend(leg_text);
        
    end
end