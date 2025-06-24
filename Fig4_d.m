clc
clear
close all

% Data definition (unit conversion: Da)
solution = {'PEO-PEG-8M', 'PEO-PEG-1M', 'PEOvis', ...
            'PS-DOP', 'PS16', 'PS7'};

Mw = [8000000, 1000000, 4000000, ...
      2000000, 16000000, 7000000] / 1000;  % Convert to kDa
b  = [40000, 2800, 20000, ...
      3000, 32400, 13700];

% Color settings
colors = [... 
    0.0, 0.4, 0.3;   % PEO-PEG-8M dark blue
    1.0, 0.5, 0.0;   % PEO-PEG-1M orange
    0.0, 0.8, 0.0;   % PEOvis       green
    0.0, 0.5, 1.0;   % PS-DOP       blue
    0.6, 0.0, 0.6;   % PS16         purple
    0.4, 0.5, 0.4];  % PS7          light green

% Plotting
figure; hold on;

for i = 1:length(Mw)
    scatter(Mw(i), b(i), 100, 'filled', ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
end

% Fitting (forced through origin)
% PEO series
idx_peo = [1 2 3];
X_peo = Mw(idx_peo)';
Y_peo = b(idx_peo)';
a_peo = (X_peo' * Y_peo) / (X_peo' * X_peo);
x_fit1 = linspace(0, max(X_peo)*1.1, 100);
y_fit1 = a_peo * x_fit1;
plot(x_fit1, y_fit1, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% PS series
idx_ps = [4 5 6];
X_ps = Mw(idx_ps)';
Y_ps = b(idx_ps)';
a_ps = (X_ps' * Y_ps) / (X_ps' * X_ps);
x_fit2 = linspace(0, max(X_ps)*1.1, 100);
y_fit2 = a_ps * x_fit2;
plot(x_fit2, y_fit2, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Legend (only scatter points, not the fit lines)
h = gobjects(1, length(solution));
for i = 1:length(solution)
    h(i) = scatter(nan, nan, 100, 'filled', ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
end
legend(h, solution, 'Location', 'southeast', 'Box', 'off', 'NumColumns', 1);

% Axis labels in LaTeX format
xlabel('$M_w$ [kDa]', 'Interpreter', 'latex');
ylabel('$b$', 'Interpreter', 'latex');

% Axis limits and style
xlim([0 20000]);
ylim([0 50000]);
set(gca, 'FontSize', 12);
% box on;
