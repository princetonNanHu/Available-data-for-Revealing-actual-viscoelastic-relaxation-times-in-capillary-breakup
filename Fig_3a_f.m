clc
clear
close all

%% ========== 1. Define fixed parameters ==========
Gamma_peo_peg_8M   = 50.26/1000;   % N/m
eta_s_peo_peg_8M   = 0.127;        % Pa·s
eta_0_peo_peg_8M   = 0.128;        % Pa·s
eta_p_peo_peg_8M   = eta_0_peo_peg_8M - eta_s_peo_peg_8M;
beta_peo_peg_8M =eta_s_peo_peg_8M/eta_0_peo_peg_8M;

R_0_exp_peo_peg = [0.975*0.470 0.500*0.49 0.263*0.3 ]/2/1000; % m
lambda_e_exp_peo_peg_drip = [5.17; 3.70; 2.53]/1000; % s
lambda_e_exp_peo_peg_drip_error = [0.29; 0.21; 0.14]/1000; % s
lambda_e_exp_peo_peg_dos = [4.93; 3.67; 2.50]/1000; % s
lambda_e_exp_peo_peg_dos_error = [0.44; 0.28; 0.40]/1000; % s

R_values  = linspace(R_0_exp_peo_peg(end)*2000,1,50) / 2 / 1000;

b_values  = [1e4, 2e4, 4e4, 8e4];

lambda_set = [5 6 7] ./ 1000;  % in seconds

%% Figure 1
% Plot layout: 3 rows, 1 column
figure('Position', [613,257,262,440]); 

%% ========== 2. Loop over three lambda values to compute and plot ==========
for iLam = 1:length(lambda_set)
    lambda_peo_peg = lambda_set(iLam);

    Ec_values = eta_p_peo_peg_8M .* R_values ./ (lambda_peo_peg .* Gamma_peo_peg_8M);
    
    results = []; % Each row: [Ec, b, min_E]

    %% ========== 2.1 Define ODE solver settings ==========
    tspan  = [0, 6000];
    RelTol = 1e-12;
    AbsTol = 1e-12;
    y0     = [1; 1; 1];   % Initial conditions: [a(0), A_zz(0), A_rr(0)]
    beta   = beta_peo_peg_8M;

    %% ========== 2.2 Loop over Ec and b combinations, solve and extract min E ==========
    for iEc = 1:length(Ec_values)
        Ec = Ec_values(iEc);

        for j = 1:length(b_values)
            b = b_values(j);

            options = odeset('Events', @stopEvent, ...
                             'RelTol', RelTol, 'AbsTol', AbsTol);
            [tSol, ySol] = ode15s(@(tau, Y) odesystem(tau, Y, Ec, beta, b), ...
                                  tspan, y0, options);

            a    = ySol(:,1);
            A_zz = ySol(:,2);
            A_rr = ySol(:,3);

            f = b ./ (b - A_zz - 2*A_rr);
            E = (1./Ec ./ a - f.*(A_zz - A_rr)) .* (1 - beta) / 3 / beta;

            min_E = Inf;
            is_decreasing = false;
            for n = 2:length(E)
                if (2/3/E(n) < 1) && (E(n) < E(n-1))
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

            results = [results; Ec, b, min_E];
        end
    end
    
    %% ========== 2.3 Plot (2/(3*minE)) vs R_0 on subplot ==========
    subplot(3, 1, iLam);
    hold on; box on;

    %% Prepare different linestyles for b values
    linestyles = {'-', '--', '-.', ':'};

    for j = 1:length(b_values)
        b_j = b_values(j);
        idx_b_j = (results(:,2) == b_j);
        Ec_j   = results(idx_b_j, 1);
        minE_j = results(idx_b_j, 3);
        [~, idx_in_EcValues] = ismember(Ec_j, Ec_values);
        R_j = R_values(idx_in_EcValues);
        [R_j_sorted, idx_sort] = sort(R_j);
        minE_j_sorted          = minE_j(idx_sort);

        exponent  = floor(log10(b_j));
        mantissa  = b_j / 10^exponent;
        label_str = sprintf('$b = %.0f \\times 10^{%d}$', mantissa, exponent);

        linestyle_idx = mod(j - 1, length(linestyles)) + 1;
        plot(R_j_sorted*2000, 2./(3.*minE_j_sorted), ...
             linestyles{linestyle_idx}, 'LineWidth', 1.5, ...
             'DisplayName', label_str,'Color','k');
    end

    x_data = R_0_exp_peo_peg;
    y_data = lambda_e_exp_peo_peg_dos/lambda_set(iLam);
    y_error = lambda_e_exp_peo_peg_dos_error/lambda_set(iLam);
    x_error = 0;
    errorbar(x_data*2000, y_data, y_error, y_error, x_error, x_error, '^', 'LineWidth', 1.2, 'CapSize', 5,'Color',[176 36 24]/256,'HandleVisibility','off');
    scatter(x_data*2000,y_data,'filled', 'MarkerFaceColor', [176 36 24]/256,'Marker','^','MarkerEdgeColor','k','HandleVisibility','off');

    x_data = R_0_exp_peo_peg;
    y_data = lambda_e_exp_peo_peg_drip/lambda_set(iLam);
    y_error = lambda_e_exp_peo_peg_drip_error/lambda_set(iLam);
    x_error = 0;
    errorbar(x_data*2000, y_data, y_error, y_error, x_error, x_error, 'diamond', 'LineWidth', 1.2, 'CapSize', 5,'Color',[176 36 24]/256,'HandleVisibility','off');
    scatter(x_data*2000,y_data,'filled', 'MarkerFaceColor', [176 36 24]/256,'Marker','diamond','MarkerEdgeColor','k','HandleVisibility','off');

    %% Replace title with inline text label
    text(0.02, 0.88, ...
         ['λ = ' num2str(lambda_peo_peg*1000) ' ms'], ...
         'Units','normalized', ...
         'FontSize',14);

    %% Show legend only in first subplot
    if iLam == 1
        legend('Location','best', 'Interpreter','latex');
    end

    %% Show x-axis label only in the third subplot
    if iLam < 3
        xlabel('');
    else
        xlabel('D_0 [mm]', 'FontSize', 14);
    end

    ylabel('λ_e / λ', 'FontSize', 14);
    set(gca, 'FontSize', 12);
end
%% Figure 3
% Plot layout: 3 rows, 1 column
figure('Position', [613,257,262,440]); 

%% ========== 2. Loop over three lambda values to compute and plot ==========
for iLam = 1:length(lambda_set)
    lambda_peo_peg = lambda_set(iLam);

    Ec_values = eta_p_peo_peg_8M .* R_values ./ (lambda_peo_peg .* Gamma_peo_peg_8M);
    
    results = []; % Each row: [Ec, b, min_E]

    %% ========== 2.1 Define ODE solver settings ==========
    tspan  = [0, 6000];
    RelTol = 1e-12;
    AbsTol = 1e-12;
    y0     = [1; 1; 1];   % Initial conditions: [a(0), A_zz(0), A_rr(0)]
    beta   = beta_peo_peg_8M;

    %% ========== 2.2 Loop over Ec and b combinations, solve and extract min E ==========
    for iEc = 1:length(Ec_values)
        Ec = Ec_values(iEc);

        for j = 1:length(b_values)
            b = b_values(j);

            options = odeset('Events', @stopEvent, ...
                             'RelTol', RelTol, 'AbsTol', AbsTol);
            [tSol, ySol] = ode15s(@(tau, Y) odesystem(tau, Y, Ec, beta, b), ...
                                  tspan, y0, options);

            a    = ySol(:,1);
            A_zz = ySol(:,2);
            A_rr = ySol(:,3);

            f = b ./ (b - A_zz - 2*A_rr);
            E = (1./Ec ./ a - f.*(A_zz - A_rr)) .* (1 - beta) / 3 / beta;

            min_E = Inf;
            is_decreasing = false;
            for n = 2:length(E)
                if (2/3/E(n) < 1) && (E(n) < E(n-1))
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

            results = [results; Ec, b, min_E];
        end
    end
    
    %% ========== 2.3 Plot normalized (2/(3*minE)) vs R_0 ==========
    subplot(3, 1, iLam);
    hold on; box on;

    %% Prepare different linestyles for b values
    linestyles = {'-', '--', '-.', ':'};

    for j = 1:length(b_values)
        b_j = b_values(j);
        idx_b_j = (results(:,2) == b_j);
        Ec_j   = results(idx_b_j, 1);
        minE_j = results(idx_b_j, 3);
        [~, idx_in_EcValues] = ismember(Ec_j, Ec_values);
        R_j = R_values(idx_in_EcValues);
        [R_j_sorted, idx_sort] = sort(R_j);
        minE_j_sorted          = minE_j(idx_sort);

        exponent  = floor(log10(b_j));
        mantissa  = b_j / 10^exponent;
        label_str = sprintf('$b = %.0f \\times 10^{%d}$', mantissa, exponent);

        linestyle_idx = mod(j - 1, length(linestyles)) + 1;
        plot(R_j_sorted*2000, 2./(3.*minE_j_sorted)/(2./(3.*minE_j_sorted(1))), ...
             linestyles{linestyle_idx}, 'LineWidth', 1.5, ...
             'DisplayName', label_str,'Color','k');
    end

    x_data = R_0_exp_peo_peg;
    y_data = lambda_e_exp_peo_peg_dos/lambda_set(iLam);
    y_error = lambda_e_exp_peo_peg_dos_error/lambda_set(iLam);
    x_error = 0;
    errorbar(x_data*2000, y_data/y_data(3), ((y_error/y_data(3)).^2+(y_data(3)*y_error(3)/y_data(3)^2)).^0.5,((y_error/y_data(3)).^2+(y_data(3)*y_error(3)/y_data(3)^2)).^0.5, x_error, x_error, '^', 'LineWidth', 1.2, 'CapSize', 5,'Color',[176 36 24]/256,'HandleVisibility','off');
    scatter(x_data*2000,y_data/y_data(3),'filled', 'MarkerFaceColor', [176 36 24]/256,'Marker','^','MarkerEdgeColor','k','HandleVisibility','off');

    x_data = R_0_exp_peo_peg;
    y_data = lambda_e_exp_peo_peg_drip/lambda_set(iLam);
    y_error = lambda_e_exp_peo_peg_drip_error/lambda_set(iLam);
    x_error = 0;
    errorbar(x_data*2000, y_data/y_data(3), ((y_error/y_data(3)).^2+(y_data(3)*y_error(3)/y_data(3)^2)).^0.5, ((y_error/y_data(3)).^2+(y_data(3)*y_error(3)/y_data(3)^2)).^0.5, x_error, x_error, 'diamond', 'LineWidth', 1.2, 'CapSize', 5,'Color',[176 36 24]/256,'HandleVisibility','off');
    scatter(x_data*2000,y_data/y_data(3),'filled', 'MarkerFaceColor', [176 36 24]/256,'Marker','diamond','MarkerEdgeColor','k','HandleVisibility','off');

    %% Replace title with inline text label
    text(0.02, 0.88, ...
         ['λ = ' num2str(lambda_peo_peg*1000) ' ms'], ...
         'Units','normalized', ...
         'FontSize',14);

    %% Show legend only in first subplot
    if iLam == 1
        legend('Location','best', 'Interpreter','latex');
    end

    %% Show x-axis label only in the third subplot
    if iLam < 3
        xlabel('');
    else
        xlabel('D_0 [mm]', 'FontSize', 14);
    end

    ylabel('λ_e / min{λ_e}', 'FontSize', 14);
    set(gca, 'FontSize', 12);
end
%% ========== Auxiliary functions: ODE system and event ==========
function dydt = odesystem(~, Y, Ec, beta, b)
    a    = Y(1);
    A_zz = Y(2);
    A_rr = Y(3);
    
    % Compute f
    f = b / (b - A_zz - 2*A_rr);
    
    % Compute E
    E = (1/Ec/a - f*(A_zz - A_rr)) * (1 - beta) / 3 / beta;
    
    % da/dtau
    da    = -0.5 * a * E;
    % dA_zz/dtau
    dA_zz =  2 * E * A_zz - (f * A_zz - 1);
    % dA_rr/dtau
    dA_rr = - E * A_rr - (f * A_rr - 1);
    
    dydt = [da; dA_zz; dA_rr];
end

function [value, isterminal, direction] = stopEvent(~, Y)
    a = Y(1);
    value      = a - 1e-9;  % Trigger event when a < 1e-9
    isterminal = 1;         % Stop integration after event
    direction  = -1;        % Detect when a crosses 1e-9 from above
end
