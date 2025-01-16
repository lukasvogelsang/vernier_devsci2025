% Code for 'The status of vernier acuity following late sight onset' 
% (Vogelsang et al., Developmental Science, 2025)

%% PART 1: SETUP
%% Folders and paths
clear all; close all; 
results_dir = 'figures';
mkdir(results_dir); 

%% Extract data
load("datafiles.mat"); % Load data
n_prakash_timepoints = size(prakash_vernier_all, 2); % Nr. of longitudinal timepoints
n_prakash = size(prakash_vernier_all, 1); % Nr. of Prakash patients

% Prakash data: mean + standard error
prakash_resacuity_mean = mean(prakash_resacuity_all, 'omitnan');
prakash_vernier_mean = mean(prakash_vernier_all, 'omitnan');
prakash_resacuity_se = std(prakash_resacuity_all, 'omitnan')./sqrt(sum(~isnan(prakash_resacuity_all)));
prakash_vernier_se = std(prakash_vernier_all, 'omitnan')./sqrt(sum(~isnan(prakash_vernier_all)));

% Control data: mean + standard error
control_200_vernier_mean = mean(control_200_vernier, 'omitnan');
control_200_vernier_se = std(control_200_vernier, 'omitnan')./sqrt(sum(~isnan(control_200_vernier)));
control_500_vernier_mean = mean(control_500_vernier, 'omitnan');
control_500_vernier_se = std(control_500_vernier, 'omitnan')./sqrt(sum(~isnan(control_500_vernier)));

%% Constants for figures
figFontSize = 12;
prakash_colors = [ 70,   6,  23;
                  157,   2,   8;
                  210,  10,   2;
                  226,  60,   4;
                  242, 130,   6;
                  255, 186,   8]./255;
control_colors = [ 0,  114, 189;
                  70,  150, 230]./255;
gray = [0.4, 0.4, 0.4];

%% PART 2: MAKE FIGURES
%% Figure 1
fig = figure('Position', [20, 70, 800, 600]); 

% Panel C (Resolution Acuity) & Panel D (Vernier Acuity)
acuity_mean = {prakash_resacuity_mean, prakash_vernier_mean};
acuity_se = {prakash_resacuity_se, prakash_vernier_se};
acuity_names = {'Resolution', 'Vernier'};
acuity_units = {'(c/deg)', '(1/min)'};
for i = 1:2
    subplot(2, 2, i);
    b = bar(1:n_prakash_timepoints, acuity_mean{i}, 'Linewidth', 1.5);
    b.FaceColor = 'flat';
    b.CData = prakash_colors(1:n_prakash_timepoints, :);
    hold on; 
    errorbar(1:n_prakash_timepoints, acuity_mean{i}, acuity_se{i}, ...
             'x', 'Color', gray, 'Linewidth', 1.5);
    title(strcat(acuity_names{i}, ' acuity means'));
    ylabel(strcat(acuity_names{i}, " acuity ", acuity_units{i}));
    xlabel('Timepoint')
    ax = gca; 
    ax.FontSize = figFontSize;
end

% Panel E: Scatter plot
subplot(2, 2, 3); 
scatter(prakash_resacuity_mean, prakash_vernier_mean, 80, ...
        prakash_colors(1:n_prakash_timepoints, :), 'filled', 'Linewidth', 2);
xlabel('Resolution acuity (c/deg)');
ylabel('Vernier acuity (1/min)'); 
title(strcat("Longitudinal means"))
ax = gca; 
ax.FontSize = figFontSize; 
grid on;
hold on;
plot([0, 7.5], [0, 0.25], '--', 'Color', gray, 'LineWidth', 2);
xlim([0, 7.5]);

% Panel F: Individual data at timepoint 6
subplot(2, 2, 4);
t = 6;
scatter(prakash_resacuity_all(:, t), prakash_vernier_all(:, t), 40, ...
        prakash_colors(t, :), 'Linewidth', 2);
title("Indiv. data at final timepoint");
xlabel('Resolution acuity (c/deg)');
ylabel('Vernier acuity (1/min)');
ax = gca;
ax.FontSize = figFontSize;
grid on;
hold on;
plot([0, 15], [0, 0.5], '--', 'Color', gray, 'LineWidth', 2);

saveas(fig, fullfile(results_dir, 'figure1.png'));

%% Supplementary Figure 1
fig = figure('Position', [20, 70, 600, 700]);
colorResacuity = [0.35, 0.35, 0.35]; 
colorVernier = [63, 152, 43]./255;

for i = 1:n_prakash
    subplot(5, 2, i);
    
    % Left y-axis 
    yyaxis left;
    plot(1:n_prakash_timepoints, prakash_resacuity_all(i, :), 'o-', ...
        'Color', colorResacuity, 'LineWidth', 2, 'MarkerFaceColor', colorResacuity);
    ylabel('Resolution');
    ylim([0 15]);
    yticks([0 5 10 15]);
    set(gca, 'ycolor', colorResacuity); 
    
    % Right y-axis 
    yyaxis right;
    plot(1:n_prakash_timepoints, prakash_vernier_all(i, :), 'd-', ...
        'Color', colorVernier, 'LineWidth', 2, 'MarkerFaceColor', colorVernier);
    ylabel('Vernier');
    ylim([0 0.5]);
    set(gca, 'ycolor', colorVernier); 
    yticks([0 0.5]);

    % Handling NaNs: plotting dashed lines between points adjacent to NaNs
    hold on;
    for k = 1:n_prakash_timepoints-1
        if isnan(prakash_resacuity_all(i, k)) && ~isnan(prakash_resacuity_all(i, k+1))
            start_idx = k;
            while isnan(prakash_resacuity_all(i, start_idx)) && start_idx > 1
                start_idx = start_idx - 1;
            end
            if start_idx < k 
                yyaxis left;
                plot([start_idx k+1], prakash_resacuity_all(i, [start_idx k+1]), ...
                     '--', 'Color', colorResacuity, 'LineWidth', 2);
            end
        end
        if isnan(prakash_vernier_all(i, k)) && ~isnan(prakash_vernier_all(i, k+1))
            start_idx = k;
            while isnan(prakash_vernier_all(i, start_idx)) && start_idx > 1
                start_idx = start_idx - 1;
            end
            if start_idx < k 
                yyaxis right;
                plot([start_idx k+1], prakash_vernier_all(i, [start_idx k+1]), ...
                     '--', 'Color', colorVernier, 'LineWidth', 2);
            end
        end
    end
    xlim([1 n_prakash_timepoints]);
    xticks(1:n_prakash_timepoints);
    if i == 9 || i == 10
        xlabel('Timepoints');
    else
        xlabel('');
    end
title(strcat("Participant ",num2str(i)))
ax = gca; ax.FontSize = 11;
end

saveas(fig, fullfile(results_dir, 'supplementary_fig1.png'));

%% Supplementary Figure 2
fig = figure('Position', [20, 70, 1250, 700]); 
x_max = ceil(max(prakash_resacuity_all, [], "all"));
y_max = ceil(max(prakash_vernier_all, [], "all")*10)/10;

for t = 1:n_prakash_timepoints
    [corrval, p_corrval] = corrcoef(prakash_vernier_all(:, t), ...
                                    prakash_resacuity_all(:, t), 'rows', 'complete');
    disp(strcat("Correlation of vernier and resolution acuity at timepoint ", ...
        num2str(t), ", p value: ", num2str(p_corrval(2))));
    subplot(2, 3, t);
    scatter(prakash_resacuity_all(:, t), prakash_vernier_all(:, t), ...
            40, prakash_colors(t, :), 'Linewidth', 2);
    hold on;
    title(strcat("Timepoint ", num2str(t), " (r = ", num2str(corrval(2)),")"))
    xlabel('Resolution acuity (c/deg)');
    ylabel('Vernier acuity (1/min)');
    ax = gca;
    ax.FontSize = figFontSize+1;
    ax.XLim = [0, x_max]; 
    ax.YLim = [0, y_max];
    grid on;
    plot([ax.XLim(1), ax.XLim(2)], [ax.XLim(1)./30, ax.XLim(2)./30], ...
         '--', 'Color', gray, 'LineWidth', 2);
end
saveas(fig, fullfile(results_dir, 'supplementary_fig2.png'));

%% Figure 2
fig = figure('Position', [20, 70, 900, 600]); 

% Panel A: Plot 20/500 comparison
subplot(2, 3, 1);
t = 3;
b = bar(1:2, [prakash_vernier_mean(t), control_500_vernier_mean], 'Linewidth', 1.5);
hold on; 
b.FaceColor = 'flat';
b.CData(1,:) = prakash_colors(t, :);
b.CData(2, :) = control_colors(1, :);
errorbar(1:2, [prakash_vernier_mean(t), control_500_vernier_mean], ...
         [prakash_vernier_se(t), control_500_vernier_se], ...
         'x', 'Color', gray, 'Linewidth', 1.5);
title('1.2 c/deg comparison');
ylabel('Vernier acuity (1/min)'); 
ax = gca;
ax.FontSize = figFontSize;
ax.YLim = [0, 0.5];
ax.XTickLabel={'Prakash','Control'};

% Panel C: Plot 20/200 comparison
subplot(2, 3, 4);
t = 6;
b = bar(1:2, [prakash_vernier_mean(t), control_200_vernier_mean], 'Linewidth', 1.5);
hold on; 
b.FaceColor = 'flat';
b.CData(1, :) = prakash_colors(t, :);
b.CData(2, :) = control_colors(2, :); 
errorbar(1:2, [prakash_vernier_mean(t), control_200_vernier_mean], ...
        [prakash_vernier_se(t), control_200_vernier_se], ...
        'x', 'Color', gray, 'Linewidth', 1.5)
title('3 c/deg comparison');
ylabel('Vernier acuity (1/min)'); 
ax = gca;
ax.FontSize = figFontSize;
ax.YLim = [0, 0.5];
ax.XTickLabel = {'Prakash','Control'};

% Panel B: Individual data (Prakash at timepoint 3 and 20/500 controls)
subplot(2, 3, 2:3);
t = 3;
non_nan_controldata = control_500_vernier(~isnan(control_500_vernier));
control_length = length(non_nan_controldata);
non_nan_prakash_data = prakash_vernier_all(~isnan(prakash_vernier_all(:, t)), t);
prakash_length = length(non_nan_prakash_data);
b = bar([non_nan_prakash_data; 0; non_nan_controldata']);
b.FaceColor = 'flat'; 
b.CData(1:prakash_length, :) = repmat(prakash_colors(t, :), prakash_length, 1); 
b.CData(prakash_length + 2:end, :) = repmat(control_colors(1, :), control_length, 1);
ax = gca; 
ax.FontSize = figFontSize; 
ax.YLim = [0, 0.5];
ylabel('Vernier acuity (1/min)'); 
title('Individual participants');

% Panel D: Individual data (Prakash at timepoint 6 and 20/200 controls)
subplot(2, 3, 5:6);
t = 6;
non_nan_controldata = control_200_vernier(~isnan(control_200_vernier));
control_length = length(non_nan_controldata);
non_nan_prakash_data = prakash_vernier_all(~isnan(prakash_vernier_all(:, t)), t);
prakash_length = length(non_nan_prakash_data);
b = bar([non_nan_prakash_data; 0; non_nan_controldata']);
b.FaceColor = 'flat'; 
b.CData(1:prakash_length, :) = repmat(prakash_colors(t, :), prakash_length, 1); 
b.CData(prakash_length + 2:end, :) = repmat(control_colors(2, :), control_length, 1);
ax = gca; 
ax.FontSize = figFontSize;
ax.YLim = [0, 0.5];
ylabel('Vernier acuity (1/min)'); 
title('Individual participants');

saveas(fig, fullfile(results_dir, 'figure2.png'));

%% Supplementary Figure 3
fig = figure('Position', [20, 70, 900, 300]); 

% Panel A: Plot 20/500 comparison (but with Prakash timepoint 2, not 3)
subplot(1, 3, 1);
t = 2;
b = bar(1:2, [prakash_vernier_mean(t), control_500_vernier_mean], 'Linewidth', 1.5);
hold on; 
b.FaceColor = 'flat';
b.CData(1,:) = prakash_colors(t, :); 
b.CData(2,:) = control_colors(1, :);
errorbar(1:2, [prakash_vernier_mean(t), control_500_vernier_mean], ...
         [prakash_vernier_se(t), control_500_vernier_se], ...
         'x', 'Color', gray, 'Linewidth', 1.5)
title('1.2 c/deg comparison');
ylabel('Vernier acuity (1/min)'); 
ax = gca;
ax.FontSize = figFontSize;
ax.YLim = [0, 0.5];
ax.XTickLabel={'Prakash','Control'};

% Panel B: Individual data (Prakash at timepoint 2 + 20/500 controls)
subplot(1, 3, 2:3);
t = 2;
non_nan_controldata = control_500_vernier(~isnan(control_500_vernier));
control_length = length(non_nan_controldata);
non_nan_prakash_data = prakash_vernier_all(~isnan(prakash_vernier_all(:, t)), t);
prakash_length = length(non_nan_prakash_data);
b = bar([non_nan_prakash_data; 0; non_nan_controldata']);
b.FaceColor = 'flat'; 
b.CData(1:prakash_length, :) = repmat(prakash_colors(t, :), prakash_length, 1); 
b.CData(prakash_length + 2:end, :) = repmat(control_colors(1, :), control_length, 1);
ax = gca; 
ax.FontSize = figFontSize; 
ax.YLim = [0, 0.5];
ylabel('Vernier acuity (1/min)'); 
title('Individual participants');

saveas(fig, fullfile(results_dir, 'supplementary_fig3.png'));

%% PART 3: STATISTICS
% Cross-group comparisons (Figure 2 + Supplementary Figure 3): 
disp(' ');
[pval] = ranksum(prakash_vernier_all(:, 2), control_500_vernier', "tail", "left"); 
disp(strcat("MWU: Prakash timepoint 2 vs. Control 20/500: p = ", num2str(pval)));
[pval] = ranksum(prakash_vernier_all(:, 3), control_500_vernier', "tail", "left"); 
disp(strcat("MWU: Prakash timepoint 3 vs. Control 20/500: p = ", num2str(pval)));
[pval] = ranksum(prakash_vernier_all(:, 6), control_200_vernier', "tail", "left");
disp(strcat("MWU: Prakash timepoint 6 vs. Control 20/200: p = ", num2str(pval)));
disp(' ');

% Stats for correlation (Figure 1E): 
[res_ver_correlation, p_res_ver_correlation] = corrcoef(prakash_resacuity_mean, prakash_vernier_mean);
disp(strcat("Correlation between the two acuity means over time: ", num2str(res_ver_correlation(2)), ...
    ", p = ", num2str(p_res_ver_correlation(2))));
disp(' ');

% Stats for Prakash vernier > resolution acuity
for timept = 1:n_prakash_timepoints
    [pval] = signrank(prakash_vernier_all(:,timept), prakash_resacuity_all(:, timept)./30, "tail", "right");
    disp(strcat("Vernier acuity > resolution acuity, for timepoint ", num2str(timept), ...
         ", p = ", num2str(pval), ", p(corrected) = ", num2str(pval*6) ));
end
disp(' ');