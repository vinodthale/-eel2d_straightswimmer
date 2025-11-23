%% ================================================================
%  NACA 4-Digit Airfoil Generator - Original/Standard Method
%
%  Implements the classic NACA 4-digit airfoil equations as described
%  in aerospace engineering textbooks (Abbott & Von Doenhoff, etc.)
%
%  Features:
%  - Standard thickness distribution (NACA 00xx series)
%  - Optional camber line (manuscript midline or custom)
%  - Perpendicular offset construction
%  - Cosine clustering for leading edge resolution
%  - High-quality visualization
%
%  Reference:
%  Abbott, I. H., & Von Doenhoff, A. E. (1959).
%  Theory of Wing Sections. Dover Publications.
% ================================================================

clear; clc; close all;

%% ---------------- User Configuration ----------------
save_plots = true;
save_coords = true;
plot_dir = fullfile(pwd, 'original_method_plots');
data_dir = fullfile(pwd, 'original_method_data');

if save_plots && ~exist(plot_dir, 'dir'), mkdir(plot_dir); end
if save_coords && ~exist(data_dir, 'dir'), mkdir(data_dir); end

%% ---------------- Airfoil Parameters ----------------
thickness_list = [0.04, 0.08, 0.12, 0.16, 0.20, 0.24];  % NACA 00xx
L = 1.0;  % chord length [m]

% NACA 4-digit thickness coefficients
a0 =  0.2969;   % √x term
a1 = -0.1260;   % x term
a2 = -0.3516;   % x² term
a3 =  0.2843;   % x³ term
a4 = -0.1015;   % x⁴ term (closed trailing edge)

%% ---------------- Camber Line Options ----------------
camber_type = 'manuscript';  % Options: 'none', 'manuscript', 'naca2412', 'custom'

% Manuscript camber line parameters
Amax = 0.125;

%% ---------------- Discretization ----------------
Nx = 256;  % number of points (cosine clustering)

% Cosine clustering (concentrates points at leading/trailing edges)
beta = linspace(0, pi, Nx);
x = 0.5 * (1 - cos(beta));  % x ∈ [0,1]

%% ---------------- Main Loop ----------------
fprintf('\n================================================================\n');
fprintf('  NACA 4-Digit Airfoil Generator - Original Method\n');
fprintf('================================================================\n\n');

for it = 1:length(thickness_list)
    t = thickness_list(it);
    tag = sprintf('%04d', round(100*t));

    fprintf('Generating NACA %s (t = %.2f)...\n', tag, t);

    %% ---- Camber Line ----
    switch camber_type
        case 'none'
            y_c = zeros(size(x));
            dyc_dx = zeros(size(x));

        case 'manuscript'
            % Manuscript midline: y_c = Amax * ((x+0.03125)/1.03125) * sin(2πx)
            y_c = Amax * ((x + 0.03125)/1.03125) .* sin(2*pi*x);
            dyc_dx = Amax * ( (1/1.03125)*sin(2*pi*x) + ...
                            ((x + 0.03125)/1.03125).*(2*pi*cos(2*pi*x)) );

        case 'naca2412'
            % Standard NACA 4-digit camber (2% max camber at 40% chord)
            p = 0.4;   % position of max camber
            m = 0.02;  % max camber

            y_c = zeros(size(x));
            dyc_dx = zeros(size(x));

            for i = 1:length(x)
                if x(i) < p
                    y_c(i) = (m/p^2) * (2*p*x(i) - x(i)^2);
                    dyc_dx(i) = (2*m/p^2) * (p - x(i));
                else
                    y_c(i) = (m/(1-p)^2) * ((1-2*p) + 2*p*x(i) - x(i)^2);
                    dyc_dx(i) = (2*m/(1-p)^2) * (p - x(i));
                end
            end

        case 'custom'
            % Add your custom camber line here
            y_c = zeros(size(x));
            dyc_dx = zeros(size(x));
    end

    %% ---- Thickness Distribution ----
    % Standard NACA formula: y_t = 5t [a0*√x + a1*x + a2*x² + a3*x³ + a4*x⁴]
    y_t = 5 * t * ( a0*sqrt(x) + a1*x + a2*x.^2 + a3*x.^3 + a4*x.^4 );
    y_t = max(y_t, 0);  % ensure non-negative

    %% ---- Surface Construction ----
    % Angle of camber line (for perpendicular offset)
    theta = atan(dyc_dx);

    % Upper surface (offset perpendicular to camber line)
    x_upper = x - y_t .* sin(theta);
    y_upper = y_c + y_t .* cos(theta);

    % Lower surface
    x_lower = x + y_t .* sin(theta);
    y_lower = y_c - y_t .* cos(theta);

    %% ---- Volume/Area Calculation ----
    % Analytical area (exact integration)
    A_analytic = trapz(x, 2*y_t);

    % Polygon area (numerical)
    px = [x_upper, fliplr(x_lower)];
    py = [y_upper, fliplr(y_lower)];
    A_polygon = abs(polyarea(px, py));

    % Error
    A_error = abs(A_polygon - A_analytic);
    A_error_pct = 100 * A_error / A_analytic;

    fprintf('  Analytical area: %.10e\n', A_analytic);
    fprintf('  Polygon area:    %.10e\n', A_polygon);
    fprintf('  Error:           %.6e (%.4f%%)\n', A_error, A_error_pct);
    fprintf('  Points:          %d\n', length(px));

    %% ---- Save Coordinates ----
    if save_coords
        % Save in standard airfoil coordinate format
        coords_file = fullfile(data_dir, sprintf('NACA%s_coords.dat', tag));
        fid = fopen(coords_file, 'w');
        fprintf(fid, 'NACA %s Airfoil\n', tag);
        fprintf(fid, '%d points\n', 2*Nx);
        fprintf(fid, 'x/c             y/c\n');

        % Standard format: start at TE upper, go to LE, then LE to TE lower
        for i = Nx:-1:1
            fprintf(fid, '%16.12f  %16.12f\n', x_upper(i), y_upper(i));
        end
        for i = 1:Nx
            fprintf(fid, '%16.12f  %16.12f\n', x_lower(i), y_lower(i));
        end
        fclose(fid);

        % Save area data
        area_file = fullfile(data_dir, sprintf('NACA%s_area.txt', tag));
        fid = fopen(area_file, 'w');
        fprintf(fid, 'NACA %s Area Analysis\n', tag);
        fprintf(fid, 'Thickness ratio: %.2f\n', t);
        fprintf(fid, 'Analytical area: %.15e\n', A_analytic);
        fprintf(fid, 'Polygon area:    %.15e\n', A_polygon);
        fprintf(fid, 'Error:           %.6e (%.4f%%)\n', A_error, A_error_pct);
        fprintf(fid, 'Number of points: %d\n', 2*Nx);
        fclose(fid);

        fprintf('  Saved: %s\n', coords_file);
    end

    %% ---- Visualization ----
    if save_plots
        fig = figure('Visible', 'off', 'Position', [100 100 1600 1000]);

        % ---- Subplot 1: Full Airfoil ----
        subplot(2,3,[1 2]);
        hold on;
        plot(x_upper, y_upper, '-', 'LineWidth', 2.5, 'Color', [0.85 0.15 0.15]);
        plot(x_lower, y_lower, '-', 'LineWidth', 2.5, 'Color', [0.15 0.15 0.85]);
        if ~strcmp(camber_type, 'none')
            plot(x, y_c, '--', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]);
        end
        plot(x_upper, y_upper, 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0.85 0.15 0.15], ...
            'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
        plot(x_lower, y_lower, 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0.15 0.15 0.85], ...
            'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
        axis equal; grid on; box on;
        xlabel('x/c', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('y/c', 'FontSize', 14, 'FontWeight', 'bold');
        title(sprintf('NACA %s Airfoil (Original Method)', tag), 'FontSize', 16, 'FontWeight', 'bold');
        if ~strcmp(camber_type, 'none')
            legend({'Upper surface', 'Lower surface', 'Camber line'}, 'Location', 'best', 'FontSize', 12);
        else
            legend({'Upper surface', 'Lower surface'}, 'Location', 'best', 'FontSize', 12);
        end
        set(gca, 'FontSize', 13);

        % ---- Subplot 2: Leading Edge Detail ----
        subplot(2,3,3);
        hold on;
        plot(x_upper, y_upper, '-o', 'LineWidth', 2.0, 'MarkerSize', 6, 'Color', [0.85 0.15 0.15]);
        plot(x_lower, y_lower, '-o', 'LineWidth', 2.0, 'MarkerSize', 6, 'Color', [0.15 0.15 0.85]);
        axis equal; grid on; box on;
        xlabel('x/c', 'FontSize', 12);
        ylabel('y/c', 'FontSize', 12);
        title('Leading Edge Detail', 'FontSize', 14, 'FontWeight', 'bold');
        xlim([-0.02 0.15]);
        ylim([-0.12 0.12]);
        set(gca, 'FontSize', 11);

        % ---- Subplot 3: Trailing Edge Detail ----
        subplot(2,3,6);
        hold on;
        plot(x_upper, y_upper, '-o', 'LineWidth', 2.0, 'MarkerSize', 6, 'Color', [0.85 0.15 0.15]);
        plot(x_lower, y_lower, '-o', 'LineWidth', 2.0, 'MarkerSize', 6, 'Color', [0.15 0.15 0.85]);
        axis equal; grid on; box on;
        xlabel('x/c', 'FontSize', 12);
        ylabel('y/c', 'FontSize', 12);
        title('Trailing Edge Detail', 'FontSize', 14, 'FontWeight', 'bold');
        xlim([0.85 1.02]);
        ylim([-0.06 0.06]);
        set(gca, 'FontSize', 11);

        % ---- Subplot 4: Thickness Distribution ----
        subplot(2,3,4);
        plot(x, 2*y_t, '-', 'LineWidth', 2.5, 'Color', [0.2 0.6 0.2]);
        grid on; box on;
        xlabel('x/c', 'FontSize', 12);
        ylabel('Thickness (2y_t/c)', 'FontSize', 12);
        title('Thickness Distribution', 'FontSize', 14, 'FontWeight', 'bold');
        set(gca, 'FontSize', 11);

        % ---- Subplot 5: Camber Line ----
        subplot(2,3,5);
        if ~strcmp(camber_type, 'none')
            plot(x, y_c, '-', 'LineWidth', 2.5, 'Color', [0.6 0.2 0.6]);
            grid on; box on;
            xlabel('x/c', 'FontSize', 12);
            ylabel('y_c/c', 'FontSize', 12);
            title('Camber Line', 'FontSize', 14, 'FontWeight', 'bold');
        else
            text(0.5, 0.5, 'No Camber\n(Symmetric Airfoil)', ...
                'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
            axis off;
        end
        set(gca, 'FontSize', 11);

        % Overall title with area info
        sgtitle(sprintf('NACA %s - Area = %.8f (Error: %.4f%%)', tag, A_polygon, A_error_pct), ...
            'FontSize', 18, 'FontWeight', 'bold');

        % Save
        outname = fullfile(plot_dir, sprintf('NACA%s_original.png', tag));
        exportgraphics(fig, outname, 'Resolution', 300);
        close(fig);

        fprintf('  Saved: %s\n', outname);
    end

    fprintf('\n');
end

fprintf('================================================================\n');
fprintf('DONE. All airfoils generated using original NACA method.\n');
if save_coords
    fprintf('Coordinate files: %s\n', data_dir);
end
if save_plots
    fprintf('Plot files:       %s\n', plot_dir);
end
fprintf('================================================================\n\n');
