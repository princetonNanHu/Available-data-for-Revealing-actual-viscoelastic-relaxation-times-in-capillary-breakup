clc
clear
close all
%% PEO_vis_PRF_2024
% give experimental Ec_e series
Ec_peo_vis_exp = (10.^[-3.278381802
-3.267022791
-3.251847983
-3.15354773
-3.122151023
-3.090546669
-3.054465997
-3.028100323
-2.981354287]);
Ec_peo_vis_exp_error =  (10.^[0.071526 0.055213 0.046789 0.047456 0.04582 0.045796 0.046244 0.045802 0.045791]);
lambda_e_peo_vis = [19.05 30.65 41.65 53.46 64.92 70.93 77.65 88.63 97.58]; %ms
lambda_e_error_peo_vis = [2.67 2.85 2.67 3.56 3.91 4.27 4.81 5.33 5.87];%ms
lambda_guess_peo_vis = 150;%ms
beta_peo_vis = 0.246/0.248;
b_e_peo_vis = [2*1e4];
%% PEG-PEO-8M
Ec_peo_peg_dos_exp = 10.^[-3.503
-3.178
-3.034];
Ec_peo_peg_dos_exp_error = 10.^([0.229
0.220
0.221]);
lambda_e_peo_peg_dos = [2.50
3.67
4.93];%ms
lambda_e_peo_peg_dos_error = [0.40
0.28
0.44];%ms
Ec_peo_peg_drip_exp = 10.^[-3.508
-3.181
-3.055];
Ec_peo_peg_drip_exp_error = 10.^[0.219
0.219
0.219];
lambda_e_peo_peg_drip = [2.53
3.70
5.17];%ms
lambda_e_peo_peg_drip_error = [0.14
0.21
0.29];%ms
lambda_guess_peo_peg = 7;%ms
b_e_peo_peg = [4*1e4];
beta_peo_peg=0.127/0.128
%% PIB-PB-0.3
Ec_pib_pb_dos_exp = 10.^[-0.995];
Ec_pib_pb_dos_exp_error = 10.^[0.047];
lambda_e_pib_pb_dos = [868];%ms
lambda_e_pib_pb_dos_error = [78.6];%ms
Ec_pib_pb_drip_exp = 10.^[-0.968 -1.298 -1.665];
Ec_pib_pb_drip_exp_error = 10.^[0.030 0.035 0.029];
lambda_e_pib_pb_drip = [815 836 831];%ms
lambda_e_pib_pb_drip_error = [32 48 28];%ms
lambda_guess_pib_pb = 900;%ms
b_e_pib_pb = [2*1e3];
beta_pib_pb=9.95/12.88
%% PS-DOP-0.05
Ec_ps_dop_dos_exp = 10.^[-1.46696253
-1.241195408
-1.131870659];
Ec_ps_dop_dos_exp_error = 10.^[0.08679417
0.053570828
0.051812343];
lambda_e_ps_dop_dos = [4.99629
4.62132
4.408596];%ms
lambda_e_ps_dop_dos_error = [0.83626
0.26507
0.21177];%ms
Ec_ps_dop_drip_exp = 10.^[-1.330253927
-1.387073395
-1.203872754
-1.152868917];
Ec_ps_dop_drip_exp_error = 10.^[0.050815009
0.058026534
0.072600976
0.049815099];
lambda_e_ps_dop_drip = [4.631146
4.156794
4.240757
4.62699];%ms
lambda_e_ps_dop_drip_error = [0.19455
0.32
0.536756
0.162359];%ms
lambda_guess_ps_dop = 4.9;%ms
b_e_ps_dop = [1e3];
beta_ps_dop=0.081/0.124
%%  PEO-PEG-1M
Ec_peo_peg_1M_dos_exp = 10.^[-1.893
-1.764
-1.546
-1.351
-1.287];
Ec_peo_peg_1M_dos_exp_error = 10.^([0.077
0.032
0.075
0.037
0.050]);
lambda_e_peo_peg_1M_dos = [1.39
1.57
1.78
1.82
2.00];%ms
lambda_e_peo_peg_1M_dos_error = [0.23
0.08
0.29
0.12
0.20];%ms

lambda_guess_peo_peg_1M = 2.2;%ms
b_e_peo_peg_1M = [7*1e2];
beta_peo_peg_1M=0.127/0.148;

%% PIB-PB-0.02

Ec_pib_pb_002_drip_exp = 10.^[-2.556
-2.538
-2.425
-2.261];
Ec_pib_pb_002_drip_exp_error = 10.^[0.108
0.088
0.091
0.097];
lambda_e_pib_pb_002_drip = [62.165
121
263.3
346.8];%ms
lambda_e_pib_pb_002_drip_error = [9.78
7.11
21.64
38.58];%ms
lambda_guess_pib_pb_002 = 700;%ms
b_e_pib_pb_002 = [2*1e3];
beta_pib_pb_002=9.95/10.22;

%% PS16-100
Ec_ps_16_100_exp = 10.^[-2.857
-2.747
-2.611];
Ec_ps_16_100_exp_error = 10.^[0.093
0.090
0.086];
lambda_e_ps_16_100= [67.97270471
79.15839537
77.13937138];%ms
lambda_e_ps_16_100_error = [6.057946678
5.961538462
2.985418211];%ms
lambda_guess_ps_16_100 = 80;%ms
b_e_ps_16_100 = 32401;
beta_ps_16_100 =0.98;
%% PS16-50
Ec_ps_16_50_exp = 10.^[-2.668
-2.523
-2.389];
Ec_ps_16_50_exp_error = 10.^[0.093
0.086
0.089];
lambda_e_ps_16_50= [43.93424318
47.17121588
46.24193548];%ms
lambda_e_ps_16_50_error = [4.044561937
2.019993751
2.981286142];%ms
lambda_guess_ps_16_50 = 50;%ms
b_e_ps_16_50 = 32401;
beta_ps_16_50 =0.98;
%% PS7-500
Ec_ps_7_500_exp = 10.^[-2.405
-2.230
-2.107];
Ec_ps_7_500_exp_error = 10.^[0.082
0.079
0.079];
lambda_e_ps_7_500= [27.97270471
28.06865178
28.16501241];%ms
lambda_e_ps_7_500_error = [2.019993751
1.443375673
1.443375673];%ms
lambda_guess_ps_7_500 = 30;%ms
b_e_ps_7_500 = 13689;
beta_ps_7_500 =0.89;
%% PS7-200
Ec_ps_7_200_exp = 10.^[-2.528
-2.432
-2.286];
Ec_ps_7_200_exp_error = 10.^[0.154
0.154
0.154];
lambda_e_ps_7_200= [15.92142266
19.15839537
18.22911497];%ms
lambda_e_ps_7_200_error = [1.159175726
1.346153846
1.251232134];%ms
lambda_guess_ps_7_200 = 20;%ms
b_e_ps_7_200 = 13689;
beta_ps_7_200 =0.96;



%% 
% Define parameter ranges (using logarithmic spacing). Ec uses the true unknown lambda
Ec_values = sort([logspace(log10(1e-6), log10(0.3), 40)]);
b_values = sort([logspace(log10(1e2), log10(1e6), 30)]); 
beta = 0.8;  % Fixed beta value

% Initialize results matrix
results = [];

% Define time span and solver settings
tspan = [0, 6000];
RelTol = 1e-12;
AbsTol = 1e-12;
a0 = 1;         % Initial value of a(0)
A_zz0 = 1;      % Initial value of A_{zz}(0)
A_rr0 = 1;      % Initial value of A_{rr}(0)
y0 = [a0; A_zz0; A_rr0];

% Loop over combinations of Ec and b
for i = 1:length(Ec_values)
    for j = 1:length(b_values)
        Ec = Ec_values(i);
        b = b_values(j);

        % Solve using ode15s
        options = odeset('Events', @stopEvent, 'RelTol', RelTol, 'AbsTol', AbsTol);
        [t, y] = ode15s(@(tau, y) odesystem(tau, y, Ec, beta, b), tspan, y0, options);

        % Extract results
        a = y(:, 1);
        A_zz = y(:, 2);
        A_rr = y(:, 3);

        % Compute E
        f = (b - 3) ./ (b - A_zz - 2 .* A_rr);
        E = (1 / Ec ./ a - f .* (A_zz - A_rr)) * (1 - beta) / 3 / beta;

        % Initialize minimum value and index
        min_E = Inf;  % Initialize with a large value
        is_decreasing = false;  % Flag to track decreasing sequence
        
        % Find the minimum value during the decreasing process
        for n = 2:length(E)
            if (2/3/E(n) < 1) && E(n) < E(n-1)
                if ~is_decreasing
                    % Start a new decreasing sequence
                    is_decreasing = true;
                    current_min_E = E(n);
                elseif E(n) < current_min_E
                    % Update the minimum in the current decreasing sequence
                    current_min_E = E(n);
                end
            elseif is_decreasing
                % Sequence has ended, check if it's the new global minimum
                is_decreasing = false;
                if current_min_E < min_E
                    min_E = current_min_E;
                end
            end
        end
        
        % Check the last decreasing sequence
        if is_decreasing && current_min_E < min_E
            min_E = current_min_E;
        end

        % Record result
        results = [results; Ec, b, min_E];
    end
end

% Convert results to grid format for 2D phase diagram plotting
[Ec_grid, b_grid] = meshgrid(Ec_values, b_values);
MinE_grid = reshape(results(:, 3), length(b_values), length(Ec_values));

% **Handle MinE_grid to avoid complex values**
MinE_grid(MinE_grid <= 0) = NaN; % Avoid log10 on zero or negative values

% Plot the 2D phase diagram
figure(1);
levels_lambda = 0:0.1:0.9;
contourf(log10(Ec_grid), log10(b_grid), (2/3 ./ MinE_grid), 100, 'Linestyle', 'none');
hold on;
contour(log10(Ec_grid), log10(b_grid), (2/3 ./ MinE_grid), levels_lambda, 'LineColor', 'k', 'Linestyle', ':');
xlabel('log_{10}(Ec=R_0\eta_p/\Gamma\lambda)');
ylabel('log_{10}(b)');
cb = colorbar;
cb.Label.String = '\lambda_{e}/\lambda';
grid on;
set(gca, 'FontSize', 14); 
cb.Label.FontSize = 14;   % Font size of colorbar label
cb.FontSize = 14;         % Font size of colorbar ticks
cb.Location = 'southoutside';  % Set colorbar position to bottom


%%
figure(4)
%%
levels_b = sort([linspace(2,5,15)]);
contourf(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), 100, 'LineColor', 'none');
hold on 
colormap(pink);
xlabel('log_{10}({Ec_e})=log_{10}(2\gamma\lambda_e/D_0\eta_p)','FontSize',14);
ylabel('log_{10}(\lambda_{e}/\lambda)','FontSize',14);
cb = colorbar;
cb.Label.String = 'log_{10}(b)';
grid on;
xlim([0.5 4]);
ylim([-1.2 0])
set(gca, 'FontSize', 14); 
cb.Label.FontSize = 14;  
cb.FontSize = 14;        
cb.Location = 'southoutside';  
%%
% PEO_vis_prf
x_data = log10(Ec_peo_vis_exp);
y_data =log10(lambda_e_peo_vis/lambda_guess_peo_vis);
y_log_error = (lambda_e_error_peo_vis./lambda_e_peo_vis)/log(10);
x_log_error = log10(Ec_peo_vis_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'o', 'LineWidth', 1.2, 'CapSize', 5,'Color',[121 121 121]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor',[121 121 121]/256,'Marker','o','MarkerEdgeColor','k');
contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_peo_vis) log10(b_e_peo_vis)], 'LineColor', [121 121 121]/256,'Linestyle','-');
%contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_peo_vis)-0.1 log10(b_e_peo_vis)+0.1], 'LineColor', [121 121 121]/256,'Linestyle','-.');
%%
% PEO_PEG_8M
x_data = log10(Ec_peo_peg_dos_exp);
y_data =log10(lambda_e_peo_peg_dos/lambda_guess_peo_peg);
y_log_error = (lambda_e_peo_peg_dos_error./lambda_e_peo_peg_dos)/log(10);
x_log_error = log10(Ec_peo_peg_dos_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, '^', 'LineWidth', 1.2, 'CapSize', 5,'Color',[176 36 24]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor', [176 36 24]/256,'Marker','^','MarkerEdgeColor','k');
x_data = log10(Ec_peo_peg_drip_exp);
y_data =log10(lambda_e_peo_peg_drip/lambda_guess_peo_peg);
y_log_error = (lambda_e_peo_peg_drip_error./lambda_e_peo_peg_drip)/log(10);
x_log_error = log10(Ec_peo_peg_drip_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'diamond', 'LineWidth', 1.2, 'CapSize', 5,'Color',[176 36 24]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor', [176 36 24]/256,'Marker','diamond','MarkerEdgeColor','k');
 contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_peo_peg) log10(b_e_peo_peg)], 'LineColor', [176 36 24]/256,'Linestyle','-');
% contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_peo_peg)-0.25 log10(b_e_peo_peg)+0.2], 'LineColor', [176 36 24]/256,'Linestyle','-.');
%%
%PIB_PB
x_data = log10(Ec_pib_pb_dos_exp);
y_data =log10(lambda_e_pib_pb_dos/lambda_guess_pib_pb);
y_log_error = (lambda_e_pib_pb_dos_error./lambda_e_pib_pb_dos)/log(10);
x_log_error = log10(Ec_pib_pb_dos_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, '^', 'LineWidth', 1.2, 'CapSize', 5,'Color',[76 123 49]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor', [76 123 49]/256,'Marker','^','MarkerEdgeColor','k');
x_data = log10(Ec_pib_pb_drip_exp);
y_data =log10(lambda_e_pib_pb_drip/lambda_guess_pib_pb);
y_log_error = (lambda_e_pib_pb_drip_error./lambda_e_pib_pb_drip)/log(10);
x_log_error = log10(Ec_pib_pb_drip_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'diamond', 'LineWidth', 1.2, 'CapSize', 5,'Color',[76 123 49]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor', [76 123 49]/256,'Marker','diamond','MarkerEdgeColor','k');
 contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_pib_pb) log10(b_e_pib_pb)], 'LineColor', [76 123 49]/256,'Linestyle','-');
% contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_pib_pb)-0.1 log10(b_e_pib_pb)+0.1], 'LineColor', [76 123 49]/256,'Linestyle','-.');

%%
%PEO_PEG_1M
x_data = log10(Ec_peo_peg_1M_dos_exp);
y_data =log10(lambda_e_peo_peg_1M_dos/lambda_guess_peo_peg_1M);
y_log_error = (lambda_e_peo_peg_1M_dos_error./lambda_e_peo_peg_1M_dos)/log(10);
x_log_error = log10(Ec_peo_peg_1M_dos_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, '^', 'LineWidth', 1.2, 'CapSize', 5,'Color',[75 41 134]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor', [75 41 134]/256,'Marker','^','MarkerEdgeColor','k');
 contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_peo_peg_1M) log10(b_e_peo_peg_1M)], 'LineColor', [75 41 134]/256,'Linestyle','-');
% contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_peo_peg_1M)-0.2 log10(b_e_peo_peg_1M)+0.2], 'LineColor', [75 41 134]/256,'Linestyle','-.');
%%
%%
%PS_DOP
x_data = log10(Ec_ps_dop_dos_exp);
y_data =log10(lambda_e_ps_dop_dos/lambda_guess_ps_dop);
y_log_error = (lambda_e_ps_dop_dos_error./lambda_e_ps_dop_dos)/log(10);
x_log_error = log10(Ec_ps_dop_dos_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, '^', 'LineWidth', 1.2, 'CapSize', 5,'Color',[34 81 150]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor', [34 81 150]/256,'Marker','^','MarkerEdgeColor','k');
x_data = log10(Ec_ps_dop_drip_exp);
y_data =log10(lambda_e_ps_dop_drip/lambda_guess_ps_dop);
y_log_error = (lambda_e_ps_dop_drip_error./lambda_e_ps_dop_drip)/log(10);
x_log_error = log10(Ec_ps_dop_drip_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'diamond', 'LineWidth', 1.2, 'CapSize', 5,'Color',[34 81 150]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor', [34 81 150]/256,'Marker','diamond','MarkerEdgeColor','k');
 contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_ps_dop) log10(b_e_ps_dop)], 'LineColor', [34 81 150]/256,'Linestyle','-');
% contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_ps_dop)-0.2 log10(b_e_ps_dop)+0.2], 'LineColor', [34 81 150]/256,'Linestyle','-.');
%%
%PIB_PB

x_data = log10(Ec_pib_pb_002_drip_exp);
y_data =log10(lambda_e_pib_pb_002_drip/lambda_guess_pib_pb_002);
y_log_error = (lambda_e_pib_pb_002_drip_error./lambda_e_pib_pb_002_drip)/log(10);
x_log_error = log10(Ec_pib_pb_002_drip_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'diamond', 'LineWidth', 1.2, 'CapSize', 5,'Color',[21 21 21]/256);
scatter(-x_data,y_data,'filled', 'MarkerFaceColor', [21 21 21]/256,'Marker','diamond','MarkerEdgeColor','k');
 contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_pib_pb_002) log10(b_e_pib_pb_002)], 'LineColor', [21 21 21]/256,'Linestyle','-');
% contour(-(log10(Ec_grid)-log10(2/3 ./ MinE_grid)),  log10(2/3 ./ MinE_grid), log10(b_grid), [log10(b_e_pib_pb_002)-0.1 log10(b_e_pib_pb_002)+0.1], 'LineColor', [21 21 21]/256,'Linestyle','-.');


%%
% PS16-100
x_data = log10(Ec_ps_16_100_exp);
y_data = log10(lambda_e_ps_16_100 / lambda_guess_ps_16_100);
y_log_error = (lambda_e_ps_16_100_error ./ lambda_e_ps_16_100) / log(10);
x_log_error = log10(Ec_ps_16_100_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'o', ...
    'LineWidth', 1.2, 'CapSize', 5, 'Color', [204 102 0]/256);
scatter(-x_data, y_data, 'filled', 'MarkerFaceColor', [204 102 0]/256, ...
    'Marker', 'o', 'MarkerEdgeColor', 'k');
contour(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), log10(2/3 ./ MinE_grid), log10(b_grid), ...
    [log10(b_e_ps_16_100) log10(b_e_ps_16_100)], 'LineColor', [204 102 0]/256, 'Linestyle', '-');
% contour(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), log10(2/3 ./ MinE_grid), log10(b_grid), ...
%     [log10(b_e_ps_16_100)-0.15 log10(b_e_ps_16_100)+0.15], 'LineColor', [204 102 0]/256, 'Linestyle', '-.');

% PS16-50
x_data = log10(Ec_ps_16_50_exp);
y_data = log10(lambda_e_ps_16_50 / lambda_guess_ps_16_50);
y_log_error = (lambda_e_ps_16_50_error ./ lambda_e_ps_16_50) / log(10);
x_log_error = log10(Ec_ps_16_50_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'o', ...
    'LineWidth', 1.2, 'CapSize', 5, 'Color', [255 178 102]/256);
scatter(-x_data, y_data, 'filled', 'MarkerFaceColor', [255 178 102]/256, ...
    'Marker', 'o', 'MarkerEdgeColor', 'k');
contour(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), log10(2/3 ./ MinE_grid), log10(b_grid), ...
     [log10(b_e_ps_16_50) log10(b_e_ps_16_50)], 'LineColor', [255 178 102]/256, 'Linestyle', '-');
% contour(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), log10(2/3 ./ MinE_grid), log10(b_grid), ...
%     [log10(b_e_ps_16_50)-0.15 log10(b_e_ps_16_50)+0.15], 'LineColor', [255 178 102]/256, 'Linestyle', '-.');

% PS7-500
x_data = log10(Ec_ps_7_500_exp);
y_data = log10(lambda_e_ps_7_500 / lambda_guess_ps_7_500);
y_log_error = (lambda_e_ps_7_500_error ./ lambda_e_ps_7_500) / log(10);
x_log_error = log10(Ec_ps_7_500_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'o', ...
    'LineWidth', 1.2, 'CapSize', 5, 'Color', [102 51 0]/256);
scatter(-x_data, y_data, 'filled', 'MarkerFaceColor', [102 51 0]/256, ...
    'Marker', 'o', 'MarkerEdgeColor', 'k');
contour(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), log10(2/3 ./ MinE_grid), log10(b_grid), ...
    [log10(b_e_ps_7_500) log10(b_e_ps_7_500)], 'LineColor', [102 51 0]/256, 'Linestyle', '-');
% contour(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), log10(2/3 ./ MinE_grid), log10(b_grid), ...
%     [log10(b_e_ps_7_500)-0.15 log10(b_e_ps_7_500)+0.15], 'LineColor', [102 51 0]/256, 'Linestyle', '-.');

% PS7-200
x_data = log10(Ec_ps_7_200_exp);
y_data = log10(lambda_e_ps_7_200 / lambda_guess_ps_7_200);
y_log_error = (lambda_e_ps_7_200_error ./ lambda_e_ps_7_200) / log(10);
x_log_error = log10(Ec_ps_7_200_exp_error);
errorbar(-x_data, y_data, y_log_error, y_log_error, x_log_error, x_log_error, 'o', ...
    'LineWidth', 1.2, 'CapSize', 5, 'Color', [153 102 51]/256);
scatter(-x_data, y_data, 'filled', 'MarkerFaceColor', [153 102 51]/256, ...
    'Marker', 'o', 'MarkerEdgeColor', 'k');
contour(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), log10(2/3 ./ MinE_grid), log10(b_grid), ...
    [log10(b_e_ps_7_200) log10(b_e_ps_7_200)], 'LineColor', [153 102 51]/256, 'Linestyle', '-');
% contour(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), log10(2/3 ./ MinE_grid), log10(b_grid), ...
%     [log10(b_e_ps_7_200)-0.15 log10(b_e_ps_7_200)+0.15], 'LineColor', [153 102 51]/256, 'Linestyle', '-.');




%%
% 定义 ODE 系统
function dydt = odesystem(~, y, Ec, beta, b)
    a = y(1);
    A_zz = y(2);
    A_rr = y(3);
    
    % 计算 f
    f = (b-3) / (b - A_zz - 2 * A_rr);
    
    % 计算 \dot{E}
    E = (1/Ec/ a - f * (A_zz - A_rr))*(1-beta)/3/beta;
    
    % 计算 da/dtau
    da = -0.5 * a * E;
    
    % 计算 dA_zz/dtau 和 dA_rr/dtau
    dA_zz = 2 * E * A_zz - (f * A_zz - 1);
    dA_rr = -E * A_rr - (f * A_rr - 1);
    
    % 将导数组合为列向量
    dydt = [da; dA_zz; dA_rr];
end

% 定义事件函数，当 a < 0 时触发停止
function [value, isterminal, direction] = stopEvent(~, y)
    value = y(1) - 1e-9;         % 检查 a 的值
    isterminal = 1;       % 当事件发生时停止积分
    direction = -1;       % 当 a 从正变负时触发事件
end