%% Reproduction of Non-Isothermal CSTR Cascades: Figures 7, 8, and 10
% Reference: Bremer & Sundmacher (2021), Front. Energy Res. 8:549298
% Equations (11)-(14) and Supplementary Equations (S22)-(S39)
clear; clc; close all;

%% 1. Parameters (Supplementary Tables S1 & S2)
p = get_reference_params();

%% ========================================================================
%  FIGURE 7: Graphical Solution of a Two-Stage CSTR Series (T_cool = 400 K)
%  Combinations [1-1], [2-1], and [3-1]
%  ========================================================================
fprintf('Generating Figure 7 (Graphical 2-stage solution at T_cool = 400 K)...\n');
n = 2;
T_cool_fig7 = 400.0; % From Fig. 7 caption: T_op^{1} = T_in = T_cool = 400 K
St_2 = p.St / n;
C_2  = p.C_1 / n;

% Find the 3 operating points (OP 1, OP 2, OP 3) of Stage 1 at T_cool = 400 K
roots_stg1 = find_stage_roots(T_cool_fig7, T_cool_fig7, 0.0, n, p);
% roots_stg1 is a struct array with fields .T and .X_stg for OP 1, 2, 3

T_grid = linspace(350, 800, 220);
X_eq_grid = arrayfun(@(T) solve_X_eq(T, p), T_grid);
% Stage 1 mass balance curve X_CO2^{1}(T)
X_mb_stg1 = arrayfun(@(T) solve_stage_mass_balance(T, 0.0, C_2, p), T_grid);
% Stage 1 energy balance line X_CO2^{1}(T)
X_eb_stg1 = ((1 + St_2) / p.dT_ad) * (T_grid - T_cool_fig7);

figure('Name', 'Figure 7: Two-Stage Graphical Solution', 'Color', 'w', 'Position', [60, 80, 1300, 430]);
titles_fig7 = {'State [1-1]', 'State [2-1]', 'State [3-1]'};

for k = 1:3
    subplot(1, 3, k); hold on; box on; set(gca, 'FontSize', 11, 'LineWidth', 1.0);
    
    % Stage 1 operating state k (OP 1, OP 2, or OP 3)
    T1_k  = roots_stg1(k).T;
    X1_k  = roots_stg1(k).X_stg;
    
    % Stage 2 parameters conditioned on Stage 1 operating point k
    dT_ad_2 = (1.0 - X1_k) * p.dT_ad;                      % Eq. (S31)
    T_op_2  = (T1_k + St_2 * T_cool_fig7) / (1.0 + St_2);  % Eq. (S19/S32)
    
    % Stage 2 mass balance curve & energy balance line
    X_mb_stg2 = arrayfun(@(T) solve_stage_mass_balance(T, X1_k, C_2, p), T_grid);
    X_eb_stg2 = ((1 + St_2) / dT_ad_2) * (T_grid - T_op_2);
    
    % Find Stage 2 OP 1 (state 1)
    roots_stg2 = find_stage_roots(T_cool_fig7, T1_k, X1_k, n, p);
    T2_1 = roots_stg2(1).T;
    X2_1 = roots_stg2(1).X_stg;
    
    % Plot Equilibrium, Stage 1 (light red/green), and Stage 2 (dark red/green)
    plot(T_grid, X_eq_grid, 'k:', 'LineWidth', 2.0);
    plot(T_grid, X_mb_stg1, 'Color', [0.96 0.45 0.42], 'LineWidth', 2.5);
    plot(T_grid, X_eb_stg1, 'Color', [0.40 0.78 0.68], 'LineWidth', 2.0);
    
    plot(T_grid, X_mb_stg2, 'Color', [0.75 0.12 0.08], 'LineWidth', 1.5);
    plot(T_grid, X_eb_stg2, 'Color', [0.05 0.45 0.40], 'LineWidth', 1.5);
    
    % Plot Stage 1 and Stage 2 Operating Point Markers
    plot(T1_k, X1_k, 'o', 'Color', [0.45 0.05 0.25], 'MarkerFaceColor', [0.75 0.65 0.75], ...
        'MarkerSize', 8, 'LineWidth', 1.5);
    plot(T2_1, X2_1, 'o', 'Color', [0.55 0.15 0.15], 'MarkerFaceColor', [0.95 0.85 0.85], ...
        'MarkerSize', 8, 'LineWidth', 1.5);
    
    text(430, 0.94, 'equilibrium', 'FontSize', 10);
    xlim([350, 800]); ylim([0, 1.08]);
    xlabel('reaction temp., T / K', 'FontSize', 11, 'FontWeight', 'bold');
    if k == 1
        ylabel('stage conversion X_{CO_2}^{\{i\}} / -', 'FontSize', 11, 'FontWeight', 'bold');
    end
    title(['Figure 7: ', titles_fig7{k}], 'FontSize', 12);
end

%% ========================================================================
%  FIGURE 8: Enumerative Solution of a Two-Stage CSTR Series (n = 2)
%  ========================================================================
fprintf('Generating Figure 8 (Two-stage bifurcation diagram, n = 2)...\n');
T_cool_vec = linspace(300, 600, 260);
paths_n2   = [1 1; 2 1; 3 1; 1 2; 2 2; 1 3; 2 3];
labels_n2  = {'1-1', '2-1', '3-1', '1-2', '2-2', '1-3', '2-3'};
colors_n2  = {[0.00 0.68 0.84], [0.55 0.02 0.25], [0.05 0.55 0.50], ...
              [0.00 0.68 0.84], [0.55 0.02 0.25], [0.05 0.55 0.50], [0.55 0.02 0.25]};
markers_n2 = {'o', 'o', 'o', '^', '^', 's', 's'};

% Store results for each of the 7 paths
X_tot_n2 = nan(size(paths_n2, 1), length(T_cool_vec));
T1_n2    = nan(size(paths_n2, 1), length(T_cool_vec));
T2_n2    = nan(size(paths_n2, 1), length(T_cool_vec));

for idx_T = 1:length(T_cool_vec)
    Tc = T_cool_vec(idx_T);
    for p_idx = 1:size(paths_n2, 1)
        state_seq = paths_n2(p_idx, :);
        [valid, X_tot, T_stg] = evaluate_cascade_path(Tc, state_seq, p);
        if valid
            X_tot_n2(p_idx, idx_T) = X_tot;
            T1_n2(p_idx, idx_T)    = T_stg(1);
            T2_n2(p_idx, idx_T)    = T_stg(2);
        end
    end
end

% Also compute the Single CSTR curve (parametric via Eq. S21) for dashed reference line
T_single_param = linspace(300, 750, 300);
X_single_param = arrayfun(@(T) solve_stage_mass_balance(T, 0.0, p.C_1, p), T_single_param);
T_op_single    = T_single_param - X_single_param * (p.dT_ad / (1.0 + p.St));

figure('Name', 'Figure 8: Two-Stage CSTR Bifurcation', 'Color', 'w', 'Position', [80, 80, 1150, 540]);

% Left Plot: Total Conversion X_CO2^{(2)} vs. T_cool
subplot(1, 2, [1]); hold on; box on; set(gca, 'FontSize', 11, 'LineWidth', 1.0);
plot(T_op_single, X_single_param, 'k--', 'LineWidth', 1.2, 'DisplayName', 'single CSTR Fig. 4');
for p_idx = 1:size(paths_n2, 1)
    plot(T_cool_vec, X_tot_n2(p_idx, :), markers_n2{p_idx}, 'Color', colors_n2{p_idx}, ...
        'MarkerSize', 4.5, 'LineWidth', 1.0, 'DisplayName', labels_n2{p_idx});
end
xlim([300, 600]); ylim([0, 1.05]);
xlabel('coolant temp. T_{cool} / K', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('total conversion X_{CO_2}^{(2)} / -', 'FontSize', 12, 'FontWeight', 'bold');
title('Figure 8 (Left): Total Conversion (n = 2)', 'FontSize', 12);
legend('Location', 'southeast', 'FontSize', 9.5);

% Top Right Plot: Stage 1 Temperature T^{1} vs. T_cool
subplot(2, 2, 2); hold on; box on; set(gca, 'FontSize', 10, 'LineWidth', 1.0);
for p_idx = 1:size(paths_n2, 1)
    plot(T_cool_vec, T1_n2(p_idx, :), markers_n2{p_idx}, 'Color', colors_n2{p_idx}, 'MarkerSize', 3.5);
end
xlim([300, 550]); ylim([300, 900]);
ylabel('reaction temp. T^{\{1\}} / K', 'FontSize', 10, 'FontWeight', 'bold');
title('Stage 1 Temperature T^{\{1\}}', 'FontSize', 11);

% Bottom Right Plot: Stage 2 Temperature T^{2} vs. T_cool
subplot(2, 2, 4); hold on; box on; set(gca, 'FontSize', 10, 'LineWidth', 1.0);
for p_idx = 1:size(paths_n2, 1)
    plot(T_cool_vec, T2_n2(p_idx, :), markers_n2{p_idx}, 'Color', colors_n2{p_idx}, 'MarkerSize', 3.5);
end
xlim([300, 550]); ylim([300, 900]);
xlabel('coolant temp. T_{cool} / K', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('reaction temp. T^{\{2\}} / K', 'FontSize', 10, 'FontWeight', 'bold');
title('Stage 2 Temperature T^{\{2\}}', 'FontSize', 11);

%% ========================================================================
%  FIGURE 10: Multi-Stage CSTR Cascades (n = 3, 4, 5, 6)
%  Uses the 4n - 1 physically valid Key-Stage Paths (Fig. 11 Hierarchy)
%  ========================================================================
fprintf('Generating Figure 10 (Multi-stage CSTR cascades n = 3, 4, 5, 6)...\n');
figure('Name', 'Figure 10: Multi-Stage CSTR Cascades (n = 3..6)', 'Color', 'w', 'Position', [70, 60, 1200, 780]);
T_cool_fig10 = linspace(200, 550, 180);
n_list = [3, 4, 5, 6];

for s_idx = 1:4
    n_stg = n_list(s_idx);
    subplot(2, 2, s_idx); hold on; box on; set(gca, 'FontSize', 10.5, 'LineWidth', 1.0);
    
    % Generate all (4*n - 1) key-stage paths from Figure 11
    [path_matrix, path_type] = generate_key_stage_paths(n_stg);
    num_paths = size(path_matrix, 1);
    X_mat = nan(num_paths, length(T_cool_fig10));
    
    for idx_T = 1:length(T_cool_fig10)
        Tc = T_cool_fig10(idx_T);
        for p_idx = 1:num_paths
            [valid, X_tot, ~] = evaluate_cascade_path(Tc, path_matrix(p_idx, :), p);
            if valid
                X_mat(p_idx, idx_T) = X_tot;
            end
        end
    end
    
    % Plot subordinate branches ("others", type == 0) in black first
    for p_idx = 1:num_paths
        if path_type(p_idx) == 0
            plot(T_cool_fig10, X_mat(p_idx, :), 'k.', 'MarkerSize', 3.5, 'HandleVisibility', 'off');
        end
    end
    plot(nan, nan, 'k.', 'MarkerSize', 6, 'DisplayName', 'others');
    
    % Plot Colored Key-Stage Inlet Branches:
    % 1: [1-1-...], 2: [2-1-...], 3: [3-1-...], 4: [2-2-...], 5: [2-3-...]
    inlet_colors  = {[0.00 0.68 0.84], [0.55 0.02 0.25], [0.05 0.55 0.50], [0.3 0.3 0.3], [0.55 0.02 0.25]};
    inlet_markers = {'.', '.', '.', 'o', '^'};
    inlet_names   = {'1-1-...', '2-1-...', '3-1-...', '2-2-...', '2-3-...'};
    for t_id = 1:5
        idx_match = find(path_type == t_id, 1);
        if strcmp(inlet_markers{t_id}, '.')
            plot(T_cool_fig10, X_mat(idx_match, :), '.', 'Color', inlet_colors{t_id}, ...
                'MarkerSize', 7, 'DisplayName', inlet_names{t_id});
        else
            plot(T_cool_fig10, X_mat(idx_match, :), inlet_markers{t_id}, 'Color', inlet_colors{t_id}, ...
                'MarkerSize', 4, 'DisplayName', inlet_names{t_id});
        end
    end
    
    xlim([200, 550]); ylim([0, 1.05]);
    xlabel('coolant temp. T_{cool} / K', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel(sprintf('total conversion X_{CO_2}^{(%d)} / -', n_stg), 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('n = %d (Total Paths = %d)', n_stg, num_paths), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'southeast', 'FontSize', 8.5);
end

%% ========================================================================
%  CORE SOLVER & HELPER FUNCTIONS
%  ========================================================================

function p = get_reference_params()
    p.R_gas      = 8.314462618;
    p.p_in       = 5.0;           % [bar]
    p.rho_cat    = 2355.2;        % [kg/m^3]
    p.eps        = 0.4;           % [-]
    p.eta_meth   = 1.0;           % [-]
    p.tau        = 0.2;           % [s]
    p.c_CO2_in   = 40.0;          % Table S1 reference [mol/m^3]
    p.St         = 3.4;           % Table S1 Stanton number [-]
    p.dT_ad      = 1074.0;        % Table S1 Adiabatic temperature rise [K]
    p.C_1        = p.tau / (p.eps * p.c_CO2_in);
    
    % Table S2 Kinetic Parameters (Koschany et al., 2016)
    p.T_ref      = 550.0;
    p.k_ref      = 3.46e-4;
    p.E_A        = 77500.0;
    p.K_OH_ref   = 0.50;          p.dH_OH  =  22400.0;
    p.K_H2_ref   = 0.44;          p.dH_H2  =  -6200.0;
    p.K_mix_ref  = 0.88;          p.dH_mix = -10000.0;
    
    p.nu         = [-1, -4, +1, +2];
    p.M          = [44.01, 2.016, 16.04, 18.015] * 1e-3;
    x_in         = [0.2, 0.8, 0.0, 0.0];
    p.w_in       = (x_in .* p.M) / sum(x_in .* p.M);
end

function [valid, X_total, T_stages] = evaluate_cascade_path(T_cool, state_seq, p)
    % Solves an n-stage cascade along a prescribed state sequence (e.g., [1, 3, 1, 1])
    n_stages = length(state_seq);
    X_cum    = 0.0;
    T_prev   = T_cool; % T_in = T_cool (Supplementary Section 2)
    T_stages = zeros(1, n_stages);
    valid    = true;
    
    for i = 1:n_stages
        req_branch = state_seq(i); % 1 = OP1 (low), 2 = OP2 (mid), 3 = OP3 (high)
        roots_i = find_stage_roots(T_cool, T_prev, X_cum, n_stages, p);
        
        if length(roots_i) < req_branch
            valid = false;
            X_total = nan;
            return;
        end
        
        T_stages(i) = roots_i(req_branch).T;
        X_cum       = X_cum + roots_i(req_branch).X_stg;
        T_prev      = T_stages(i);
    end
    X_total = X_cum;
end

function roots_struct = find_stage_roots(T_cool, T_prev, X_cum_in, n_stages, p)
    % Finds all steady-state roots (OP 1, OP 2, OP 3) of Stage {i}
    St_n    = p.St / n_stages;
    C_n     = p.C_1 / n_stages;
    dT_ad_i = (1.0 - X_cum_in) * p.dT_ad;             % Eq. (S31)
    T_op_i  = (T_prev + St_n * T_cool) / (1.0 + St_n);% Eq. (S19/S32)
    
    if dT_ad_i <= 1e-7
        roots_struct = struct('T', T_op_i, 'X_stg', 0.0);
        return;
    end
    
    % Maximum possible stage conversion cannot exceed (1 - X_cum_in)
    T_max = T_op_i + dT_ad_i / (1.0 + St_n);
    
    % Stage residual in terms of T_i: X_eb(T_i) - C_n * R_meth(T_i, X_cum_in + X_eb(T_i)) = 0
    res_func = @(T) stage_residual(T, T_op_i, dT_ad_i, St_n, X_cum_in, C_n, p);
    
    % Coarse scan across [T_op_i, T_max] to bracket all sign changes (up to 3 roots)
    N_scan = 85;
    T_scan = linspace(T_op_i + 1e-7, T_max - 1e-6, N_scan);
    f_scan = arrayfun(res_func, T_scan);
    
    roots_struct = struct('T', {}, 'X_stg', {});
    r_count = 0;
    for k = 1:(N_scan - 1)
        if f_scan(k) * f_scan(k+1) <= 0
            T_root = fzero(res_func, [T_scan(k), T_scan(k+1)]);
            X_root = ((1.0 + St_n) / dT_ad_i) * (T_root - T_op_i);
            r_count = r_count + 1;
            roots_struct(r_count).T     = T_root;
            roots_struct(r_count).X_stg = max(0.0, X_root);
        end
    end
end

function res = stage_residual(T, T_op_i, dT_ad_i, St_n, X_cum_in, C_n, p)
    X_stg = ((1.0 + St_n) / dT_ad_i) * (T - T_op_i); % Eq. (S32)
    X_tot = X_cum_in + X_stg;
    
    % Beyond chemical equilibrium, R_meth is negative
    R_val = R_meth_signed(X_tot, T, p);
    res   = X_stg - C_n * R_val;
end

function [path_matrix, path_type] = generate_key_stage_paths(n)
    % Generates the (4n - 1) physically significant paths from Figure 11
    path_matrix = zeros(4*n - 1, n);
    path_type   = zeros(4*n - 1, 1);
    
    % 1. All extinguished: [1-1-...-1]
    path_matrix(1, :) = ones(1, n);
    path_type(1)      = 1; % Type 1: 1-1-...
    idx = 1;
    
    % 2. Single stage in State 2: [...-1-2-1-...] (n paths)
    for i = 1:n
        idx = idx + 1;
        row = ones(1, n); row(i) = 2;
        path_matrix(idx, :) = row;
        if i == 1, path_type(idx) = 2; end % Type 2: 2-1-...
    end
    
    % 3. Single stage in State 3: [...-1-3-1-...] (n paths)
    for i = 1:n
        idx = idx + 1;
        row = ones(1, n); row(i) = 3;
        path_matrix(idx, :) = row;
        if i == 1, path_type(idx) = 3; end % Type 3: 3-1-...
    end
    
    % 4. Two adjacent stages in [2-2]: [...-1-2-2-1-...] (n-1 paths)
    for i = 1:(n - 1)
        idx = idx + 1;
        row = ones(1, n); row(i) = 2; row(i+1) = 2;
        path_matrix(idx, :) = row;
        if i == 1, path_type(idx) = 4; end % Type 4: 2-2-...
    end
    
    % 5. Two adjacent stages in [2-3]: [...-1-2-3-1-...] (n-1 paths)
    for i = 1:(n - 1)
        idx = idx + 1;
        row = ones(1, n); row(i) = 2; row(i+1) = 3;
        path_matrix(idx, :) = row;
        if i == 1, path_type(idx) = 5; end % Type 5: 2-3-...
    end
end

function X_stg = solve_stage_mass_balance(T, X_cum_in, C_n, p)
    X_eq = solve_X_eq(T, p);
    if X_cum_in >= X_eq - 1e-11
        X_stg = 0.0;
        return;
    end
    res = @(X_s) X_s - C_n * R_meth_signed(X_cum_in + X_s, T, p);
    X_stg = fzero(res, [0.0, X_eq - X_cum_in]);
end

function X_eq = solve_X_eq(T, p)
    X_eq = fzero(@(X) R_meth_signed(X, T, p), [1e-12, 1.0 - 1e-12]);
end

function R_meth = R_meth_signed(X_CO2_cum, T, p)
    X_CO2_cum = max(0.0, min(X_CO2_cum, 1.0 - 1e-14));
    w = -X_CO2_cum .* (p.nu ./ p.nu(1)) .* (p.M ./ p.M(1)) .* p.w_in(1) + p.w_in;
    moles_rel = w ./ p.M;
    x_alpha   = moles_rel ./ sum(moles_rel);
    p_alpha   = x_alpha .* p.p_in;
    
    p_CO2 = p_alpha(1); p_H2 = p_alpha(2); p_CH4 = p_alpha(3); p_H2O = p_alpha(4);
    
    inv_T = (1.0 / p.T_ref) - (1.0 / T);
    k     = p.k_ref     * exp((p.E_A    / p.R_gas) * inv_T);
    K_OH  = p.K_OH_ref  * exp((p.dH_OH  / p.R_gas) * inv_T);
    K_H2  = p.K_H2_ref  * exp((p.dH_H2  / p.R_gas) * inv_T);
    K_mix = p.K_mix_ref * exp((p.dH_mix / p.R_gas) * inv_T);
    K_eq  = 137.0 * (T ^ (-3.998)) * exp(158700.0 / (p.R_gas * T));
    
    DEN    = 1.0 + K_OH * (p_H2O / sqrt(p_H2)) + K_H2 * sqrt(p_H2) + K_mix * sqrt(p_CO2);
    r_meth = k * sqrt(p_CO2 * p_H2) * (1.0 - (p_CH4 * p_H2O^2) / (K_eq * p_CO2 * p_H2^4)) / (DEN^2);
    R_meth = (1.0 - p.eps) * p.eta_meth * (r_meth * p.rho_cat * 1000.0);
end