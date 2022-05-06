function FIGURE_NozzlePerformance
% This function plots some specific results for the performance of
% homogenised models for activation travelling through a nozzle-type
% structure

% Specify the type of homogenisation for which results are plotted
homog_type = 'linear';

% Specify the results to plot (varying only r (1D) or l and r (2D))
homog_extents1D = [5, 10, 25, 50];
homog_extent2D = 10;

% Specify the left mouth width to use for the 1D plot
plot_l = 200;

% Specify colours
colors1D = [ [0.1 0.1 0.7];
    [0.45 0.45 0.95];
    [0.95 0.45 0.45];
    [0.8 0.15 0.15]
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the data for the no-fibrosis case and calculate the activation
% time of the end opposite stimulus in this case
load('funnel_nofib.mat','nofib_act_map');
nofib_AT = nanmean(nofib_act_map(:,end));

% Load in the finescale data and convert widths to actual measurements
load('funnel_finescale.mat','fine_act_maps','widths_LR');
widths_LR = 10*widths_LR;

% Load in the homogenised data and process it
load(['funnelRG_',homog_type,'.mat'], 'homog_act_maps', 'homog_extents');

% Determine which homogenisation level the given extent corresponds to
homog_level = find( homog_extent2D == homog_extents );
if isempty(homog_level)
    error('Homogenisation extent given does not appear to match the results collected');
end




% FIGURE 6B: RESULTS VARYING ONE PARAMETER (R)

% A dummy value used for when the simulation blocks
block_delay = 2.6;
block_jitter = 0.1;

% Create a new figure
figure('units','normalized','OuterPosition',[0 0 1 1]);  hold on;

% Caculate the delay for each case
finescale_AD = zeros(size(widths_LR,1),1);
for m = 1:size(widths_LR,1)
    finescale_AT = nanmean(fine_act_maps{m}(:,end));
    if finescale_AT == -1
        finescale_AD(m) = NaN;
    else
        finescale_AD(m) = finescale_AT - nofib_AT;
    end
end

% Plot the delay for this case
nanvals = (isnan(finescale_AD) & widths_LR(:,1) == plot_l);
nonnanvals = (~isnan(finescale_AD) & widths_LR(:,1) == plot_l);
plot( widths_LR(nonnanvals,2), finescale_AD(nonnanvals), 'k.', 'MarkerSize', 50);
plot( widths_LR(nanvals,2), block_delay, 'kx', 'MarkerSize', 20, 'LineWidth', 2.5);

for n = 1:size(homog_act_maps,1)
    
    % Caculate the delay for each case
    homog_AD = zeros(size(widths_LR,1),1);
    for m = 1:size(widths_LR,1)
        homog_AT = nanmean(homog_act_maps{n,m}(:,end));
        if homog_AT == -1
            homog_AD(m) = NaN;
        else
            homog_AD(m) = homog_AT - nofib_AT;
        end
    end
    
    nanvals = (isnan(homog_AD) & widths_LR(:,1) == plot_l);
    nonnanvals_normal = (~isnan(homog_AD) & widths_LR(:,1) == plot_l & homog_AD <= 3);
    nonnanvals_extreme = (~isnan(homog_AD) & widths_LR(:,1) == plot_l & homog_AD > 3);
    plot( widths_LR(nonnanvals_normal,2), homog_AD(nonnanvals_normal), '.', 'MarkerSize', 40, 'MarkerEdgeColor', colors1D(n,:));
    if any(nonnanvals_extreme)
        plot( widths_LR(nonnanvals_extreme,2), block_delay + n*block_jitter, 'O', 'MarkerSize', 15, 'MarkerEdgeColor', colors1D(n,:), 'LineWidth', 2.5);
    end
    plot( widths_LR(nanvals,2), block_delay + n*block_jitter, 'x', 'MarkerSize', 15, 'MarkerEdgeColor', colors1D(n,:), 'LineWidth', 2.5);
    
end

% Plot a dotted line to deliniate block versus delay
plot([-10, 510],[2.5,2.5],'--','LineWidth',2,'Color',[0.4 0.4 0.4])
xlim([0 500]);

% Set up the x axis
xlabel('Exit Mouth Width, r (\mum)');

% Set up the y axis - custom ticks to mark block
ylabel('Activation Delay (ms)');
ylim([0 3.1]);
yticks([0 1 2 2.8]);
yticklabels({'0','1','2','BLOCK'});

% Set fontsize
set(gca,'FontSize', 24);

% Set square axis
axis square;

% Add legend
leg_txt = {'Finescale (10 \mum)'};
leg_pos = [7, 12,18,26];
for k = 1:length(findobj(gca,'type','line'))
    if ismember(k,leg_pos)
        leg_txt{k} = ['Homogenised (',num2str(homog_extents1D(k==leg_pos)*10),' \mum)'];
    else
        leg_txt{k} = '';
    end
end
leg_obj = legend(leg_txt,'Location','northeastoutside');
pos = get(leg_obj,'Position');
set(leg_obj,'Position',[pos(1)*0.925, pos(2)-0.05, pos(3), pos(4)]);




% FIGURE 6C: RESULTS VARYING TWO PARAMETERS (L AND R)

% Create a figure
figure('units','normalized','OuterPosition',[0 0 1 1]);

% Finescale
for m = 1:size(widths_LR,1)
    
    % Activation time at end of domain
    finescale_AT = nanmean(fine_act_maps{m}(:,end));
    % Convert to a delay measure
    if finescale_AT == -1
        finescale_AD(m) = NaN;
    else
        finescale_AD(m) = finescale_AT - nofib_AT;
    end
    % Activation time at end of domain
    homog_AT = nanmean(homog_act_maps{homog_level,m}(:,end));
    % Convert to a delay measure
    if homog_AT == -1
        homog_AD(m) = NaN;
    else
        homog_AD(m) = homog_AT - nofib_AT;
    end
    
end

% Create a colourmap for showing the extent of delay
load('extra_colormaps.mat', 'viridis');
lr_clr = viridis;

% Prepare figure
subplot(1,2,1); hold on;
% Plot delay values as colours
scatter( widths_LR(:,1), widths_LR(:,2), 100, finescale_AD, 'filled' );
% Plot X's for cases of block
plot( widths_LR(isnan(finescale_AD),1), widths_LR(isnan(finescale_AD),2), 'x', 'MarkerEdgeColor',[0 0 0], 'MarkerSize', 20, 'LineWidth', 3);
% Specify colour axis and apply colormap
caxis([0 2]);
colormap(lr_clr);
% Specify axis limits
axis([0 520 0 520]);
% Add axis labels and set font
xlabel('Entrance Mouth Width, l (\mu m)');
ylabel('Exit Mouth Width, r (\mu m)');
set(gca,'FontSize',24);
title('Finescale (10 \mum)', 'FontSize', 28);

% Prepare figure
subplot(1,2,2); hold on;
% Plot delay values as colours
scatter( widths_LR(:,1), widths_LR(:,2), 100, homog_AD, 'filled' );
% Plot X's for cases of block
plot( widths_LR(isnan(homog_AD),1), widths_LR(isnan(homog_AD),2), 'x', 'MarkerEdgeColor',[0 0 0], 'MarkerSize', 20, 'LineWidth', 3);
% Specify colour axis
caxis([0 2]);
colormap(lr_clr);
% Specify axis limits
axis([0 520 0 520]);
% Add colorbar
colorbar_obj=colorbar;
colorbar_obj.Position = colorbar_obj.Position + [0.05 0 0 0];
title(colorbar_obj, 'AD (ms)', 'FontSize', 24);
% Add axis labels and set font
xlabel('Entrance Mouth Width, l (\mu m)');
ylabel('Exit Mouth Width, r (\mu m)');
set(gca,'FontSize',24);
title(['Homogenised (',num2str(10*homog_extent2D),' \mum)'],'FontSize', 28);



%%% FIGURE 6D: EXAMPLE ACTIVATION TIME MAPS

% Specify the problems that are to be plotted
problems = [10, 102];

% Create a colormap using 'plasma' with appended colours for obstruction
% and non-activated sites
load('extra_colormaps.mat','plasma');
act_clr = [ [0.1 0.1 0.1]; plasma; [0.6 0.6 0.6] ];

% Initialise figure
figure('units', 'normalized', 'OuterPosition',[0 0 1 0.6]);

% Create the axes
ax_properties.xgap = 0.025;
ax_properties.margin = 0.05;
ax_properties.ygap = 0.04;
ax_properties.leftspace = 0.1;
ax_properties.bottomspace = 0.1;
ax_properties.topspace = 0.1;
ax_properties.rightspace = 0.1;
ax = createAxes(4, 2, 2, ax_properties);

% Loop over the list of problems, plotting each in turn.
global_max = 0;
for k = 1:length(problems)
   
    % Read out the activation maps
    fine_AT = fine_act_maps{problems(k)};
    homog_AT = homog_act_maps{homog_level,problems(k)};
      
    % Use an extreme negative value to mark the fibrosis
    fine_AT( isnan(fine_AT) ) = -100;
    homog_AT( isnan(homog_AT) ) = -100;
    % Use an extreme positive value to mark sites that failed to activate
    fine_AT( fine_AT == -1 ) = 100;
    homog_AT( homog_AT == -1 ) = 100; 
    
    % Create a list of 'valid' sites (i.e. those that are not indicating fibrosis or block)
    fine_valid = ( fine_AT(:) ~= -100 & fine_AT(:) ~= 100 );
    homog_valid = ( homog_AT(:) ~= -100 & homog_AT(:) ~= 100 );
    
    % Subtract off the time before activation at the start (using
    % finescale)
    start_val = nanmean(fine_AT(:,1));
    fine_AT(fine_valid) = fine_AT(fine_valid) - start_val;
    homog_AT(homog_valid) = homog_AT(homog_valid) - start_val;
    
    % Find limits of plot (consistent colour axis)
    maxval = max( [fine_AT(fine_valid); homog_AT(homog_valid)] );
    global_max = max([ maxval, global_max ]);
        
    % Plot these
    imagesc(ax{1,k},fine_AT);
    axis(ax{1,k},'equal');
    colormap(ax{1,k},act_clr);
    
    imagesc(ax{2,k},homog_AT);
    axis(ax{2,k},'equal');
    colormap(ax{2,k},act_clr);
    
end

% Loop over all plots and set them to a single consistent colour axis
for k = 1:length(problems)
    
    %subplot(length(problems),2,sub2ind([length(problems),2],k,1));
    caxis(ax{1,k},[-global_max*0.1, global_max*1.1]);
    
    % Add labels if leftmost map
    if k == 1
       yticks(ax{1,k},size(fine_AT,1) / 2);
       yticklabels(ax{1,k},{'Finescale   '});
       set(ax{1,k},'FontSize', 24);
    else
        yticks(ax{1,k},[]);
        yticklabels(ax{1,k},[]);
    end
    xticks(ax{1,k},[]);
    xticklabels(ax{1,k},[]);
    xlim(ax{1,k},[0.5 size(fine_AT,2)-0.5]);
    ylim(ax{1,k},[0.5 size(fine_AT,1)-0.5]);
    
    %subplot(length(problems),2,sub2ind([length(problems),2],k,2));
    caxis(ax{2,k},[-global_max*0.1, global_max*1.1]);
    
    % Add labels if leftmost map
    if k == 1
       yticks(ax{2,k},size(homog_AT,1) / 2);
       yticklabels(ax{2,k},{'Homogenised   '});
       set(ax{2,k},'FontSize', 24);
    else
        yticks(ax{2,k},[]);
        yticklabels(ax{2,k},[]);
    end
    xticks(ax{2,k},[]);
    xticklabels(ax{2,k},[]);
    xlim(ax{2,k},[0.5 size(homog_AT,2)-0.5]);
    ylim(ax{2,k},[0.5 size(homog_AT,1)-0.5]);
    
end

% Add highlight circles onto the first homogenised model
rectangle(ax{2,1},'Position',[0.485*size(homog_AT,2), 0.26*size(homog_AT,1), 0.225*size(homog_AT,1), 0.225*size(homog_AT,1)],'Curvature',[1 1], 'EdgeColor', [0 1 1], 'LineWidth', 3);
rectangle(ax{2,1},'Position',[0.485*size(homog_AT,2), 0.60*size(homog_AT,1), 0.225*size(homog_AT,1), 0.225*size(homog_AT,1)],'Curvature',[1 1], 'EdgeColor', [0 1 1], 'LineWidth', 3);

% Add a colorbar
margin = 0.05;
titleSpace = 0.05;
axpos = get(ax{2,2},'Position');
clrbar_obj = colorbar(ax{2,2},'Position', [axpos(1)+axpos(3)*17.5/16,  axpos(2)+margin, 0.015, 1-8*margin - titleSpace], 'FontSize', 24);
title(clrbar_obj,'AT (ms)');
