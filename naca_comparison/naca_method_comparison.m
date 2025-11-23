%% ================================================================
%  NACA Airfoil Method Comparison
%  Compares:
%      METHOD 1: Original NACA 4-digit formulation (Standard Textbook)
%      METHOD 2: IBAMR Perpendicular Stacking Method (Current Implementation)
%
%  Reference: Standard NACA airfoil equations
%  - NACA 4-digit series thickness distribution
%  - Symmetric airfoils (00xx series)
%  - Perpendicular offset from camber line
%
%  Output:
%  - Side-by-side comparison plots
%  - Detailed geometric analysis
%  - Volume/area verification
%  - Point distribution analysis
% ================================================================

clear; clc; close all;

%% ---------------- Configuration ----------------
save_plots = true;
% Use centralized output directory
output_base = fullfile(fileparts(pwd), 'output');
plot_dir = fullfile(output_base, 'images');
if save_plots && ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

%% ---------------- Thickness Set ----------------
thickness_list = [0.04, 0.08, 0.12, 0.16, 0.20, 0.24];

%% ---------------- NACA 4-digit Coefficients ----------------
% Standard NACA thickness distribution coefficients
a0 =  0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 =  0.2843;
a4 = -0.1015;  % closed trailing edge

%% ---------------- Chord and Discretization ----------------
L = 1.0;  % chord length [m]

% High-resolution for original method (continuous surface)
Nx_original = 256;

% IBAMR stacking parameters (matching manuscript)
domain_x = 10.0 * L;
domain_y = 4.0 * L;
coarse_patch_x = 40;
coarse_patch_y = 16;
MAX_LEVELS = 4;
REF_RATIO = 4;

dx_coarse = domain_x / coarse_patch_x;
dy_coarse = domain_y / coarse_patch_y;
dx_fine = dx_coarse / (REF_RATIO^(MAX_LEVELS-1));
dy_fine = dy_coarse / (REF_RATIO^(MAX_LEVELS-1));
ds = dx_fine;  % layer spacing for stacking

Nx_stack = max(128, ceil(L / dx_fine));

%% ---------------- Camber Line (Manuscript) ----------------
Amax = 0.125;
use_midline = true;

%% ---------------- Placement Parameters ----------------
Xcenter = domain_x/2 + 0.48;
Ycenter = domain_y/2;

%% ---------------- Comparison Table Storage ----------------
comparison_table = [];

fprintf('\n');
fprintf('================================================================\n');
fprintf('  NACA Airfoil Method Comparison\n');
fprintf('================================================================\n\n');
fprintf('METHOD 1: Original NACA (Continuous Surface)\n');
fprintf('  - Standard textbook formulation\n');
fprintf('  - Perpendicular thickness distribution from camber line\n');
fprintf('  - Nx = %d points\n\n', Nx_original);
fprintf('METHOD 2: IBAMR Stacking (Discrete Layers)\n');
fprintf('  - Perpendicular stacking with ds = %.6e\n', ds);
fprintf('  - IBAMR-style layering (j=0..k)\n');
fprintf('  - Nx = %d chordwise sections\n\n', Nx_stack);
fprintf('================================================================\n\n');

%% ================================================================
%  MAIN COMPARISON LOOP
%% ================================================================
for it = 1:length(thickness_list)
    t = thickness_list(it);
    tag = sprintf('%04d', round(100*t));

    fprintf('Processing NACA %s (t = %.2f)...\n', tag, t);

    %% ============================================================
    %  METHOD 1: Original NACA Formulation
    %% ============================================================

    % Chordwise parametric coordinate (non-uniform, clustered at LE)
    beta = linspace(0, pi, Nx_original);
    x_orig = 0.5 * (1 - cos(beta));  % cosine clustering

    % Camber line
    if use_midline
        y_camber_orig = Amax * ((x_orig + 0.03125)/1.03125) .* sin(2*pi*x_orig);
        dyc_dx_orig = Amax * ( (1/1.03125)*sin(2*pi*x_orig) + ...
                              ((x_orig + 0.03125)/1.03125).*(2*pi*cos(2*pi*x_orig)) );
    else
        y_camber_orig = zeros(size(x_orig));
        dyc_dx_orig = zeros(size(x_orig));
    end

    % Thickness distribution
    yt_orig = 5 * t * ( a0*sqrt(x_orig) + a1*x_orig + a2*x_orig.^2 + ...
                        a3*x_orig.^3 + a4*x_orig.^4 );
    yt_orig = max(yt_orig, 0);

    % Perpendicular offset angle
    theta_orig = atan(dyc_dx_orig);

    % Upper and lower surfaces (perpendicular to camber line)
    xu_orig = x_orig - yt_orig .* sin(theta_orig);
    yu_orig = y_camber_orig + yt_orig .* cos(theta_orig);

    xl_orig = x_orig + yt_orig .* sin(theta_orig);
    yl_orig = y_camber_orig - yt_orig .* cos(theta_orig);

    % Analytical area (exact NACA formula)
    Vf_analytic_orig = trapz(x_orig, 2.0 * yt_orig);

    % Surface polygon area
    px_orig = [xu_orig, fliplr(xl_orig)];
    py_orig = [yu_orig, fliplr(yl_orig)];
    Vf_surface_orig = abs(polyarea(px_orig, py_orig));

    % Number of surface points
    Npts_orig = 2 * Nx_original;

    %% ============================================================
    %  METHOD 2: IBAMR Perpendicular Stacking
    %% ============================================================

    % Uniform chordwise distribution
    x_stack = linspace(0, L, Nx_stack);

    % Storage for all stacked vertices
    Coord_stack = cell(Nx_stack, 2);
    Ntot_stack = 0;

    % Storage for outermost envelope
    xu_stack_env = zeros(1, Nx_stack);
    yu_stack_env = zeros(1, Nx_stack);
    xl_stack_env = zeros(1, Nx_stack);
    yl_stack_env = zeros(1, Nx_stack);

    % Storage for analytical surface on same grid
    xu_stack_surf = zeros(1, Nx_stack);
    yu_stack_surf = zeros(1, Nx_stack);
    xl_stack_surf = zeros(1, Nx_stack);
    yl_stack_surf = zeros(1, Nx_stack);

    max_k = 0;

    for i = 1:Nx_stack
        x_phys = x_stack(i);
        x_norm = min(1.0, x_phys / L);

        % Camber line
        if use_midline
            y_c = Amax * ((x_phys + 0.03125)/1.03125) * sin(2*pi*x_phys);
            dyc_dx = Amax * ( (1/1.03125)*sin(2*pi*x_phys) + ...
                            ((x_phys + 0.03125)/1.03125)*(2*pi*cos(2*pi*x_phys)) );
        else
            y_c = 0.0;
            dyc_dx = 0.0;
        end

        % Thickness
        yt = 5.0 * t * ( a0*sqrt(max(x_norm,0)) + a1*x_norm + ...
                         a2*x_norm^2 + a3*x_norm^3 + a4*x_norm^4 );
        yt = max(yt, 0.0);

        % Normal direction (perpendicular to camber)
        theta = atan(dyc_dx);
        nx = -sin(theta);
        ny = cos(theta);

        % Analytical surface points (for comparison)
        base = [x_phys - 0.5; y_c];
        upp_surf = base + yt * [nx; ny];
        low_surf = base - yt * [nx; ny];
        upp_world = upp_surf + [Xcenter; Ycenter];
        low_world = low_surf + [Xcenter; Ycenter];
        xu_stack_surf(i) = upp_world(1);
        yu_stack_surf(i) = upp_world(2);
        xl_stack_surf(i) = low_world(1);
        yl_stack_surf(i) = low_world(2);

        % Determine number of layers
        k = floor(yt / ds);
        if k > max_k, max_k = k; end

        NumUp = max(1, k+1);
        NumDown = k;

        x_up = zeros(1, NumUp);
        y_up = zeros(1, NumUp);
        x_down = zeros(1, NumDown);
        y_down = zeros(1, NumDown);

        % Stack upper layers (j=0..k)
        for j = 0:(NumUp-1)
            offset = j * ds;
            pt = base + offset * [nx; ny];
            pt_world = pt + [Xcenter; Ycenter];
            x_up(j+1) = pt_world(1);
            y_up(j+1) = pt_world(2);
            Ntot_stack = Ntot_stack + 1;
        end

        % Stack lower layers (j=1..k)
        for j = 1:NumDown
            offset = j * ds;
            pt = base - offset * [nx; ny];
            pt_world = pt + [Xcenter; Ycenter];
            x_down(j) = pt_world(1);
            y_down(j) = pt_world(2);
            Ntot_stack = Ntot_stack + 1;
        end

        % Store all vertices
        Coord_stack{i,1} = [x_up, x_down];
        Coord_stack{i,2} = [y_up, y_down];

        % Outermost envelope
        outer_up = base + k * [nx; ny];
        outer_low = base - k * [nx; ny];
        outer_up_world = outer_up + [Xcenter; Ycenter];
        outer_low_world = outer_low + [Xcenter; Ycenter];
        xu_stack_env(i) = outer_up_world(1);
        yu_stack_env(i) = outer_up_world(2);
        xl_stack_env(i) = outer_low_world(1);
        yl_stack_env(i) = outer_low_world(2);
    end

    % Surface area from analytical surface on stacking grid
    px_stack_surf = [xu_stack_surf, fliplr(xl_stack_surf)];
    py_stack_surf = [yu_stack_surf, fliplr(yl_stack_surf)];
    Vf_surface_stack = abs(polyarea(px_stack_surf, py_stack_surf));

    % Envelope area
    px_stack_env = [xu_stack_env, fliplr(xl_stack_env)];
    py_stack_env = [yu_stack_env, fliplr(yl_stack_env)];
    Vf_envelope_stack = abs(polyarea(px_stack_env, py_stack_env));

    %% ============================================================
    %  COMPUTE DIFFERENCES
    %% ============================================================

    % Area differences
    delta_Vf_orig = abs(Vf_surface_orig - Vf_analytic_orig);
    delta_Vf_stack = abs(Vf_surface_stack - Vf_analytic_orig);
    delta_Vf_env = abs(Vf_envelope_stack - Vf_analytic_orig);

    % Point count ratio
    point_ratio = Ntot_stack / Npts_orig;

    %% ============================================================
    %  STORE COMPARISON DATA
    %% ============================================================

    comparison_table = [comparison_table; struct(...
        't', t, ...
        'Vf_analytic', Vf_analytic_orig, ...
        'Vf_orig_surface', Vf_surface_orig, ...
        'Vf_stack_surface', Vf_surface_stack, ...
        'Vf_stack_envelope', Vf_envelope_stack, ...
        'delta_orig', delta_Vf_orig, ...
        'delta_stack_surf', delta_Vf_stack, ...
        'delta_stack_env', delta_Vf_env, ...
        'Npts_orig', Npts_orig, ...
        'Npts_stack', Ntot_stack, ...
        'point_ratio', point_ratio, ...
        'max_k', max_k)];

    fprintf('  Original Method: %d points, Vf = %.9e, error = %.3e\n', ...
        Npts_orig, Vf_surface_orig, delta_Vf_orig);
    fprintf('  Stacking Method: %d points, Vf = %.9e, error = %.3e\n', ...
        Ntot_stack, Vf_surface_stack, delta_Vf_stack);
    fprintf('  Point ratio: %.2f, max layers: %d\n\n', point_ratio, max_k);

    %% ============================================================
    %  VISUALIZATION
    %% ============================================================

    if save_plots
        fig = figure('Visible', 'off', 'Position', [100 100 1800 1200]);

        % ---- Subplot 1: Original Method (Full Domain) ----
        subplot(2,3,1);
        hold on;
        plot(xu_orig, yu_orig, '-', 'LineWidth', 2.5, 'Color', [0.85 0.15 0.15]);
        plot(xl_orig, yl_orig, '-', 'LineWidth', 2.5, 'Color', [0.85 0.15 0.15]);
        if use_midline
            plot(x_orig, y_camber_orig, '--', 'LineWidth', 1.5, 'Color', [0.4 0.4 0.4]);
        end
        axis equal; grid on; box on;
        xlabel('x/c', 'FontSize', 12);
        ylabel('y/c', 'FontSize', 12);
        title(sprintf('METHOD 1: Original NACA %s', tag), 'FontSize', 14, 'FontWeight', 'bold');
        legend({'Upper surface', 'Lower surface', 'Camber line'}, 'Location', 'best', 'FontSize', 10);
        set(gca, 'FontSize', 11);

        % ---- Subplot 2: Stacking Method (Full Domain) ----
        subplot(2,3,2);
        hold on;
        % Plot all stacked points
        for i = 1:Nx_stack
            plot(Coord_stack{i,1}, Coord_stack{i,2}, '.', 'MarkerSize', 4, 'Color', [0.1 0.1 0.1]);
        end
        % Analytical surface overlay
        plot(xu_stack_surf, yu_stack_surf, '-', 'LineWidth', 2.0, 'Color', [0.85 0.15 0.15]);
        plot(xl_stack_surf, yl_stack_surf, '-', 'LineWidth', 2.0, 'Color', [0.85 0.15 0.15]);
        % Envelope
        plot(xu_stack_env, yu_stack_env, '--', 'LineWidth', 1.5, 'Color', [0.1 0.4 0.8]);
        plot(xl_stack_env, yl_stack_env, '--', 'LineWidth', 1.5, 'Color', [0.1 0.4 0.8]);
        axis equal; grid on; box on;
        xlabel('x [m]', 'FontSize', 12);
        ylabel('y [m]', 'FontSize', 12);
        title(sprintf('METHOD 2: IBAMR Stacking NACA %s', tag), 'FontSize', 14, 'FontWeight', 'bold');
        legend({'Stacked points', 'Analytical surface', '', 'Outer envelope'}, 'Location', 'best', 'FontSize', 10);
        xlim([Xcenter-0.6 Xcenter+1.2]);
        ylim([Ycenter-0.4 Ycenter+0.4]);
        set(gca, 'FontSize', 11);

        % ---- Subplot 3: Leading Edge Detail (Original) ----
        subplot(2,3,4);
        hold on;
        plot(xu_orig, yu_orig, '-', 'LineWidth', 2.5, 'Color', [0.85 0.15 0.15]);
        plot(xl_orig, yl_orig, '-', 'LineWidth', 2.5, 'Color', [0.85 0.15 0.15]);
        plot(xu_orig, yu_orig, 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.85 0.15 0.15], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        plot(xl_orig, yl_orig, 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.85 0.15 0.15], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        axis equal; grid on; box on;
        xlabel('x/c', 'FontSize', 12);
        ylabel('y/c', 'FontSize', 12);
        title('Leading Edge Detail (Original)', 'FontSize', 13, 'FontWeight', 'bold');
        xlim([-0.05 0.2]);
        ylim([-0.15 0.15]);
        set(gca, 'FontSize', 11);

        % ---- Subplot 4: Leading Edge Detail (Stacking) ----
        subplot(2,3,5);
        hold on;
        for i = 1:Nx_stack
            plot(Coord_stack{i,1}, Coord_stack{i,2}, '.', 'MarkerSize', 6, 'Color', [0.1 0.1 0.1]);
        end
        plot(xu_stack_surf, yu_stack_surf, '-', 'LineWidth', 2.0, 'Color', [0.85 0.15 0.15]);
        plot(xl_stack_surf, yl_stack_surf, '-', 'LineWidth', 2.0, 'Color', [0.85 0.15 0.15]);
        axis equal; grid on; box on;
        xlabel('x [m]', 'FontSize', 12);
        ylabel('y [m]', 'FontSize', 12);
        title('Leading Edge Detail (Stacking)', 'FontSize', 13, 'FontWeight', 'bold');
        xlim([Xcenter-0.55 Xcenter-0.35]);
        ylim([Ycenter-0.1 Ycenter+0.1]);
        set(gca, 'FontSize', 11);

        % ---- Subplot 5: Thickness Distribution Comparison ----
        subplot(2,3,3);
        hold on;
        plot(x_orig, 2*yt_orig, '-', 'LineWidth', 2.5, 'Color', [0.85 0.15 0.15]);
        yt_stack_vec = zeros(1, Nx_stack);
        for i = 1:Nx_stack
            x_norm = min(1.0, x_stack(i) / L);
            yt_stack_vec(i) = 5.0 * t * ( a0*sqrt(max(x_norm,0)) + a1*x_norm + ...
                                          a2*x_norm^2 + a3*x_norm^3 + a4*x_norm^4 );
        end
        plot(x_stack, 2*yt_stack_vec, 'o', 'MarkerSize', 6, 'LineWidth', 1.5, 'Color', [0.1 0.4 0.8]);
        grid on; box on;
        xlabel('x/c', 'FontSize', 12);
        ylabel('Thickness (2y_t)', 'FontSize', 12);
        title('Thickness Distribution', 'FontSize', 13, 'FontWeight', 'bold');
        legend({'Original (continuous)', 'Stacking (discrete)'}, 'Location', 'best', 'FontSize', 10);
        set(gca, 'FontSize', 11);

        % ---- Subplot 6: Area Comparison ----
        subplot(2,3,6);
        hold on;
        bar_data = [Vf_analytic_orig, Vf_surface_orig, Vf_surface_stack, Vf_envelope_stack];
        bar(bar_data, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
        set(gca, 'XTickLabel', {'Analytic', 'Orig Surf', 'Stack Surf', 'Stack Env'});
        ylabel('Area [m^2]', 'FontSize', 12);
        title('Volume/Area Comparison', 'FontSize', 13, 'FontWeight', 'bold');
        grid on; box on;
        xtickangle(45);
        set(gca, 'FontSize', 11);

        % Add text annotations
        for i = 1:length(bar_data)
            text(i, bar_data(i)*1.02, sprintf('%.6f', bar_data(i)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        end

        % Overall title
        sgtitle(sprintf('NACA %s Airfoil Method Comparison (t = %.2f)', tag, t), ...
            'FontSize', 16, 'FontWeight', 'bold');

        % Save figure
        outname = fullfile(plot_dir, sprintf('comparison_NACA%s.png', tag));
        exportgraphics(fig, outname, 'Resolution', 300);
        close(fig);
        fprintf('  Saved: %s\n\n', outname);
    end
end

%% ================================================================
%  PRINT COMPREHENSIVE COMPARISON TABLE
%% ================================================================

fprintf('\n================================================================\n');
fprintf('  COMPREHENSIVE COMPARISON TABLE\n');
fprintf('================================================================\n\n');

fprintf('Part 1: Volume/Area Analysis\n');
fprintf('-----------------------------------------------------------------------------------\n');
fprintf('  t     Vf_analytic      Vf_orig         Vf_stack_surf   Vf_stack_env\n');
fprintf('-----------------------------------------------------------------------------------\n');
for i = 1:length(comparison_table)
    r = comparison_table(i);
    fprintf('%5.2f   %13.9e   %13.9e   %13.9e   %13.9e\n', ...
        r.t, r.Vf_analytic, r.Vf_orig_surface, r.Vf_stack_surface, r.Vf_stack_envelope);
end
fprintf('-----------------------------------------------------------------------------------\n\n');

fprintf('Part 2: Error Analysis\n');
fprintf('-------------------------------------------------------------------------\n');
fprintf('  t     Error_orig       Error_stack     Error_env       Error_ratio\n');
fprintf('-------------------------------------------------------------------------\n');
for i = 1:length(comparison_table)
    r = comparison_table(i);
    err_ratio = r.delta_stack_surf / (r.delta_orig + eps);
    fprintf('%5.2f   %13.6e   %13.6e   %13.6e   %10.3f\n', ...
        r.t, r.delta_orig, r.delta_stack_surf, r.delta_stack_env, err_ratio);
end
fprintf('-------------------------------------------------------------------------\n\n');

fprintf('Part 3: Computational Efficiency\n');
fprintf('---------------------------------------------------------------\n');
fprintf('  t     Npts_orig   Npts_stack   Point_ratio   Max_layers\n');
fprintf('---------------------------------------------------------------\n');
for i = 1:length(comparison_table)
    r = comparison_table(i);
    fprintf('%5.2f   %8d   %10d   %11.2f   %10d\n', ...
        r.t, r.Npts_orig, r.Npts_stack, r.point_ratio, r.max_k);
end
fprintf('---------------------------------------------------------------\n\n');

%% ================================================================
%  KEY FINDINGS SUMMARY
%% ================================================================

fprintf('================================================================\n');
fprintf('  KEY FINDINGS\n');
fprintf('================================================================\n\n');

avg_error_orig = mean([comparison_table.delta_orig]);
avg_error_stack = mean([comparison_table.delta_stack_surf]);
avg_point_ratio = mean([comparison_table.point_ratio]);

fprintf('1. ACCURACY:\n');
fprintf('   - Original method avg error:  %.6e\n', avg_error_orig);
fprintf('   - Stacking method avg error:  %.6e\n', avg_error_stack);
fprintf('   - Relative accuracy: %.2f%%\n\n', (avg_error_stack/avg_error_orig)*100);

fprintf('2. COMPUTATIONAL COST:\n');
fprintf('   - Average point ratio (stack/orig): %.2f\n', avg_point_ratio);
fprintf('   - Stacking uses ~%.0f%% more points on average\n\n', (avg_point_ratio-1)*100);

fprintf('3. METHOD CHARACTERISTICS:\n');
fprintf('   Original Method:\n');
fprintf('     + Smooth continuous surface representation\n');
fprintf('     + Direct implementation of NACA equations\n');
fprintf('     + Minimal points for surface definition\n');
fprintf('     - No internal structure for fluid coupling\n\n');
fprintf('   Stacking Method:\n');
fprintf('     + Full 3D volume representation (IBAMR compatible)\n');
fprintf('     + Internal layers for fluid-structure interaction\n');
fprintf('     + Mesh-adaptive (ds = dx_fine)\n');
fprintf('     - Higher point count\n');
fprintf('     - Discretization-dependent accuracy\n\n');

fprintf('================================================================\n');
fprintf('DONE. Comparison complete.\n');
fprintf('Plots saved to: %s\n', plot_dir);
fprintf('================================================================\n\n');
