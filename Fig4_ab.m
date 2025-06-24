clc;
clear;
close all;

%% PEO_vis_PRF_2024
% Provide experimental Ec_e series
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
beta_peo_peg = 0.127/0.128;

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
beta_pib_pb = 9.95/12.88;

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
beta_ps_dop = 0.081/0.124;

%% PEO-PEG-1M
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
beta_peo_peg_1M = 0.127/0.148;
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
beta_pib_pb_002 = 9.95/10.22;

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
beta_ps_16_100 = 0.98;

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
beta_ps_16_50 = 0.98;

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
beta_ps_7_500 = 0.89;

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
beta_ps_7_200 = 0.96;

%% Define parameter ranges (using logarithmic spacing)
Ec_values = sort([logspace(log10(1e-6), log10(0.3), 20)]);
b_values = sort([logspace(log10(1e2), log10(1e6), 20)]); 

% Beta values to iterate over (linear scale)
beta_values = [beta_peo_vis, beta_peo_peg, beta_pib_pb, beta_ps_dop, beta_peo_peg_1M, beta_pib_pb_002];
beta_values = [beta_values, beta_ps_16_100, beta_ps_16_50, beta_ps_7_500, beta_ps_7_200];

% Generate mesh grid
[Ec_grid, b_grid] = meshgrid(Ec_values, b_values);

% Initialize solver parameters
tspan = [0, 6000];
RelTol = 1e-12;
AbsTol = 1e-12;
a0 = 1;         
A_zz0 = 1;      
A_rr0 = 1;      
y0 = [a0; A_zz0; A_rr0];

% Create 3D surface figure
figure(3);
hold on;
colormap("copper");
view(3);
for k = 1:length(beta_values)
    beta = beta_values(k);
    results = [];
    
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
    
            % Calculate strain rate E
            f = (b - 3) ./ (b - A_zz - 2 .* A_rr);
            E = (1 / Ec ./ a - f .* (A_zz - A_rr)) * (1 - beta) / (3 * beta);
    
            % Find the minimum of E when it's decreasing and 2/3E < 1
            min_E = Inf;
            is_decreasing = false;
            
            for n = 2:length(E)
                if (2/3/E(n) < 1) && E(n) < E(n-1)
                    if ~is_decreasing
                        is_decreasing = true;
                        current_min_E = E(n);
                    elseif E(n) < current_min_E
                        current_min_E = E(n);
                    end
                elseif is_decreasing
                    is_decreasing = false;
                    if current_min_E < min_E
                        min_E = current_min_E;
                    end
                end
            end
    
            if is_decreasing && current_min_E < min_E
                min_E = current_min_E;
            end

            % Store the result
            results = [results; Ec, b, min_E];
        end
    end

    % Convert results to grid format
    MinE_grid = reshape(results(:, 3), length(b_values), length(Ec_values));

    % Handle non-positive values to avoid complex numbers
    MinE_grid(MinE_grid <= 0) = NaN; % Avoid log10 of negative or zero

    % Draw 3D surface plot
    surf(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), ...
         log10(2/3 ./ MinE_grid), ...
         ones(size(log10(Ec_grid))) * beta, ...
         log10(b_grid), ...
         'EdgeColor', 'none');
    alpha(0.1);  % Set transparency
end

xlabel('log_{10}({Ec_e})=log_{10}(2\gamma\lambda_e/D_0\eta_p)', 'FontSize', 16);
ylabel('log_{10}(\lambda_{e}/\lambda)', 'FontSize', 16);
zlabel('\beta', 'FontSize', 16);  % Beta is not log-transformed
cb = colorbar;
cb.Label.String = 'log_{10}(b)';
grid on;
set(gca, 'FontSize', 14); 
cb.Label.FontSize = 16;
cb.FontSize = 16;

% Set axis display range
xlim([0.5 4]);          % Range for log(Ec_e)
ylim([-1.2, 0]);        % Range for log(lambda_e/lambda)
hold off;
hold on;

% Iterate through different beta layers
for k = 1:length(beta_values)
    beta = beta_values(k);
    
    % Handle different experimental datasets
    if k == 1
        % PEO_vis_prf
        x_data_dos = log10(Ec_peo_vis_exp);
        y_data_dos = log10(lambda_e_peo_vis / lambda_guess_peo_vis);
        y_log_error_dos = (lambda_e_error_peo_vis ./ lambda_e_peo_vis) / log(10);
        x_log_error_dos = log10(Ec_peo_vis_exp_error);
        color = [121 121 121] / 256;
        b_value = log10(b_e_peo_vis);
        
        % No drip data
        x_data_drip = [];
        y_data_drip = [];
    elseif k == 2
        % PEO_PEG_8M
        x_data_dos = log10(Ec_peo_peg_dos_exp);
        y_data_dos = log10(lambda_e_peo_peg_dos / lambda_guess_peo_peg);
        y_log_error_dos = (lambda_e_peo_peg_dos_error ./ lambda_e_peo_peg_dos) / log(10);
        x_log_error_dos = log10(Ec_peo_peg_dos_exp_error);

        x_data_drip = log10(Ec_peo_peg_drip_exp);
        y_data_drip = log10(lambda_e_peo_peg_drip / lambda_guess_peo_peg);
        y_log_error_drip = (lambda_e_peo_peg_drip_error ./ lambda_e_peo_peg_drip) / log(10);
        x_log_error_drip = log10(Ec_peo_peg_drip_exp_error);

        color = [176 36 24] / 256;
        b_value = log10(b_e_peo_peg);
    elseif k == 3
        % PIB_PB
        x_data_dos = log10(Ec_pib_pb_dos_exp);
        y_data_dos = log10(lambda_e_pib_pb_dos / lambda_guess_pib_pb);
        y_log_error_dos = (lambda_e_pib_pb_dos_error ./ lambda_e_pib_pb_dos) / log(10);
        x_log_error_dos = log10(Ec_pib_pb_dos_exp_error);

        x_data_drip = log10(Ec_pib_pb_drip_exp);
        y_data_drip = log10(lambda_e_pib_pb_drip / lambda_guess_pib_pb);
        y_log_error_drip = (lambda_e_pib_pb_drip_error ./ lambda_e_pib_pb_drip) / log(10);
        x_log_error_drip = log10(Ec_pib_pb_drip_exp_error);

        color = [76 123 49] / 256;
        b_value = log10(b_e_pib_pb);
    elseif k == 4
        % PEO_PEG_1M
        x_data_dos = log10(Ec_peo_peg_1M_dos_exp);
        y_data_dos = log10(lambda_e_peo_peg_1M_dos / lambda_guess_peo_peg_1M);
        y_log_error_dos = (lambda_e_peo_peg_1M_dos_error ./ lambda_e_peo_peg_1M_dos) / log(10);
        x_log_error_dos = log10(Ec_peo_peg_1M_dos_exp_error);

        % No drip data
        x_data_drip = [];
        y_data_drip = [];

        color = [75 41 134] / 256;
        b_value = log10(b_e_peo_peg_1M);
    elseif k == 5
        % PS_DOP
        x_data_dos = log10(Ec_ps_dop_dos_exp);
        y_data_dos = log10(lambda_e_ps_dop_dos / lambda_guess_ps_dop);
        y_log_error_dos = (lambda_e_ps_dop_dos_error ./ lambda_e_ps_dop_dos) / log(10);
        x_log_error_dos = log10(Ec_ps_dop_dos_exp_error);

        x_data_drip = log10(Ec_ps_dop_drip_exp);
        y_data_drip = log10(lambda_e_ps_dop_drip / lambda_guess_ps_dop);
        y_log_error_drip = (lambda_e_ps_dop_drip_error ./ lambda_e_ps_dop_drip) / log(10);
        x_log_error_drip = log10(Ec_ps_dop_drip_exp_error);

        color = [34 81 150] / 256;
        b_value = log10(b_e_ps_dop);
    elseif k == 6
        % PIB_PB_002
        x_data_dos = [];
        y_data_dos = [];

        x_data_drip = log10(Ec_pib_pb_002_drip_exp);
        y_data_drip = log10(lambda_e_pib_pb_002_drip / lambda_guess_pib_pb_002);
        y_log_error_drip = (lambda_e_pib_pb_002_drip_error ./ lambda_e_pib_pb_002_drip) / log(10);
        x_log_error_drip = log10(Ec_pib_pb_002_drip_exp_error);

        color = [21 21 21] / 256;
        b_value = log10(b_e_pib_pb_002);
    elseif k == 7
        % PS_16_100
        x_data_dos = log10(Ec_ps_16_100_exp);
        y_data_dos = log10(lambda_e_ps_16_100 / lambda_guess_ps_16_100);
        y_log_error_dos = (lambda_e_ps_16_100_error ./ lambda_e_ps_16_100) / log(10);
        x_log_error_dos = log10(Ec_ps_16_100_exp_error);

        x_data_drip = [];
        y_data_drip = [];

        color = [204 102 0] / 256;  % Dark orange
        b_value = log10(b_e_ps_16_100);
    elseif k == 8
        % PS_16_50
        x_data_dos = log10(Ec_ps_16_50_exp);
        y_data_dos = log10(lambda_e_ps_16_50 / lambda_guess_ps_16_50);
        y_log_error_dos = (lambda_e_ps_16_50_error ./ lambda_e_ps_16_50) / log(10);
        x_log_error_dos = log10(Ec_ps_16_50_exp_error);

        x_data_drip = [];
        y_data_drip = [];

        color = [255 178 102] / 256;  % Light orange
        b_value = log10(b_e_ps_16_50);
    elseif k == 9
        % PS_7_500
        x_data_dos = log10(Ec_ps_7_500_exp);
        y_data_dos = log10(lambda_e_ps_7_500 / lambda_guess_ps_7_500);
        y_log_error_dos = (lambda_e_ps_7_500_error ./ lambda_e_ps_7_500) / log(10);
        x_log_error_dos = log10(Ec_ps_7_500_exp_error);

        x_data_drip = [];
        y_data_drip = [];

        color = [102 51 0] / 256;  % Dark brown
        b_value = log10(b_e_ps_7_500);
    elseif k == 10
        % PS_7_200
        x_data_dos = log10(Ec_ps_7_200_exp);
        y_data_dos = log10(lambda_e_ps_7_200 / lambda_guess_ps_7_200);
        y_log_error_dos = (lambda_e_ps_7_200_error ./ lambda_e_ps_7_200) / log(10);
        x_log_error_dos = log10(Ec_ps_7_200_exp_error);

        x_data_drip = [];
        y_data_drip = [];

        color = [153 102 51] / 256;  % Light brown
        b_value = log10(b_e_ps_7_200);
    end
    % Plot DOS data (square symbols)
    scatter3(-x_data_dos, y_data_dos, ones(size(x_data_dos)) * beta, ...
             'filled', 'MarkerFaceColor', color, 'Marker', '^', 'MarkerEdgeColor', 'k');

    % Plot DRIP data (diamond symbols)
    scatter3(-x_data_drip, y_data_drip, ones(size(x_data_drip)) * beta, ...
             'filled', 'MarkerFaceColor', color, 'Marker', 'diamond', 'MarkerEdgeColor', 'k');

    % Plot error bars for DOS data
    for i = 1:length(x_data_dos)
        % X error bars (Ec)
        line([-x_data_dos(i) - x_log_error_dos(i), -x_data_dos(i) + x_log_error_dos(i)], ...
             [y_data_dos(i), y_data_dos(i)], ...
             [beta, beta], 'Color', color, 'LineWidth', 1.2);
        % Y error bars (lambda_e)
        line([-x_data_dos(i), -x_data_dos(i)], ...
             [y_data_dos(i) - y_log_error_dos(i), y_data_dos(i) + y_log_error_dos(i)], ...
             [beta, beta], 'Color', color, 'LineWidth', 1.2);
    end

    % Plot error bars for DRIP data
    for i = 1:length(x_data_drip)
        % X error bars (Ec)
        line([-x_data_drip(i) - x_log_error_drip(i), -x_data_drip(i) + x_log_error_drip(i)], ...
             [y_data_drip(i), y_data_drip(i)], ...
             [beta, beta], 'Color', color, 'LineWidth', 1.2);
        % Y error bars (lambda_e)
        line([-x_data_drip(i), -x_data_drip(i)], ...
             [y_data_drip(i) - y_log_error_drip(i), y_data_drip(i) + y_log_error_drip(i)], ...
             [beta, beta], 'Color', color, 'LineWidth', 1.2);
    end

    % Plot contour line corresponding to b_value (solid)
    contour3(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), ...
             log10(2/3 ./ MinE_grid), ...
             ones(size(log10(Ec_grid))) * beta, ...
             log10(b_grid), [b_value b_value], ...
             'LineColor', color, 'LineStyle', '-', 'LineWidth', 1.5);  

    % Plot b ± 0.1 contour (dash-dot)
    contour3(-(log10(Ec_grid) - log10(2/3 ./ MinE_grid)), ...
             log10(2/3 ./ MinE_grid), ...
             ones(size(log10(Ec_grid))) * beta, ...
             log10(b_grid), [b_value - 0.1, b_value + 0.1], ...
             'LineColor', color, 'LineStyle', '-.', 'LineWidth', 1.5);
end

hold off;
colorbar('off');
%% Define the ODE system
function dydt = odesystem(~, y, Ec, beta, b)
    a = y(1);
    A_zz = y(2);
    A_rr = y(3);
    
    % Compute f
    f = (b - 3) / (b - A_zz - 2 * A_rr);
    
    % Compute strain rate E
    E = (1 / Ec / a - f * (A_zz - A_rr)) * (1 - beta) / (3 * beta);
    
    % Compute da/dtau
    da = -0.5 * a * E;
    
    % Compute dA_zz/dtau and dA_rr/dtau
    dA_zz = 2 * E * A_zz - (f * A_zz - 1);
    dA_rr = -E * A_rr - (f * A_rr - 1);
    
    % Return the derivative vector
    dydt = [da; dA_zz; dA_rr];
end

%% Define the event function: stop integration when a < 0
function [value, isterminal, direction] = stopEvent(~, y)
    value = y(1) - 1e-9;       % Monitor when a gets very small
    isterminal = 1;            % Stop integration when event occurs
    direction = -1;            % Trigger only when a is decreasing
end
