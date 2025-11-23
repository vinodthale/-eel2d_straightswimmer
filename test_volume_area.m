%% ================================================================
%  Volume and Area Verification Test for Eel2D Geometry
%
%  Comprehensive validation of geometric properties including:
%  1. Total body area (cross-sectional in 2D)
%  2. Volume calculation (2D area × unit depth)
%  3. Section-wise area distribution (head vs tail)
%  4. Cross-sectional height verification
%  5. Comparison with analytical formulas
%  6. Conservation checks
%% ================================================================

clear; clc; close all;

fprintf('\n');
fprintf('================================================================\n');
fprintf('  EEL2D VOLUME AND AREA VERIFICATION TEST\n');
fprintf('================================================================\n\n');

%% ---------------- Configuration (matching main script) ----------------
% Background Mesh Properties
Ny = 16*4*4*4; Ly = 4;
Nx = 32*4*4*4; Lx = 8;
dy = Ly/Ny; dx = Lx/Nx;

L  = 1;                 % length of the body
Wh = 0.04*L;            % width of the head
Xh = 0.04;

% No of points on backbone
BodyNx = ceil(L/dx);
HeadNx = ceil(0.04/dx);
Coord = cell(BodyNx, 2);
NumMatPoints = 0;

angleRotation  = 0;
RotationMatrix = [cos(angleRotation) -sin(angleRotation);
                  sin(angleRotation)  cos(angleRotation)];

fprintf('Configuration:\n');
fprintf('  Body length (L): %.2f\n', L);
fprintf('  Head width (Wh): %.4f\n', Wh);
fprintf('  Head length (Xh): %.4f\n', Xh);
fprintf('  Grid spacing: dx = %.6e, dy = %.6e\n', dx, dy);
fprintf('  Sections: HeadNx = %d, BodyNx = %d\n\n', HeadNx, BodyNx);

%% ================================================================
%  PART 1: GENERATE GEOMETRY AND STORE ENVELOPE
%% ================================================================
fprintf('Generating eel geometry...\n');

% Arrays to store envelope (upper and lower surfaces)
x_envelope_upper = zeros(1, BodyNx);
y_envelope_upper = zeros(1, BodyNx);
x_envelope_lower = zeros(1, BodyNx);
y_envelope_lower = zeros(1, BodyNx);
height_actual = zeros(1, BodyNx);
height_theoretical = zeros(1, BodyNx);
area_cross_section = zeros(1, BodyNx);

% Store from nose to head of the body
for i = 1:HeadNx
    x = (i-1)*dx;
    y = 0.125 * ((x + 0.03125)/1.03125) * sin(2*pi*x);
    height = sqrt(2*Wh*x - x^2);
    NumPointsInHeight = ceil(height/dy);

    ycoord_up = []; ycoord_down = []; xcoord_up = []; xcoord_down = [];
    for j = 1:NumPointsInHeight
        xshifted = x - 0.5;
        yshifted = y;
        RotatedCoord = RotationMatrix*[xshifted; yshifted];

        xcoord_up(j)   = RotatedCoord(1) - (j-1)*dy*sin(angleRotation);
        xcoord_down(j) = RotatedCoord(1) + j*dy*sin(angleRotation);
        ycoord_up(j)   = RotatedCoord(2) + (j-1)*dy*cos(angleRotation);
        ycoord_down(j) = RotatedCoord(2) - j*dy*cos(angleRotation);
        NumMatPoints = NumMatPoints + 2;
    end

    Coord{i,1} = cat(2, xcoord_up, xcoord_down);
    Coord{i,2} = cat(2, ycoord_up, ycoord_down);

    % Store envelope
    x_envelope_upper(i) = xcoord_up(end);
    y_envelope_upper(i) = ycoord_up(end);
    x_envelope_lower(i) = xcoord_down(end);
    y_envelope_lower(i) = ycoord_down(end);

    % Calculate actual height from discretization
    if NumPointsInHeight > 0
        height_actual(i) = (NumPointsInHeight - 1) * dy;
    else
        height_actual(i) = 0;
    end
    height_theoretical(i) = height;

    % Cross-sectional area (2 * half-height for symmetric shape)
    area_cross_section(i) = 2 * height;
end

% Store from head to tail of the body
for i = HeadNx+1:BodyNx
    x = (i-1)*dx;
    y = 0.125 * ((x + 0.03125)/1.03125) * sin(2*pi*x);
    height = Wh * (L - x) / (L - Xh);
    NumPointsInHeight = ceil(height/dy);

    ycoord_up = []; ycoord_down = []; xcoord_up = []; xcoord_down = [];
    for j = 1:NumPointsInHeight
        xshifted = x - 0.5;
        yshifted = y;
        RotatedCoord = RotationMatrix*[xshifted; yshifted];

        xcoord_up(j)   = RotatedCoord(1) - (j-1)*dy*sin(angleRotation);
        xcoord_down(j) = RotatedCoord(1) + j*dy*sin(angleRotation);
        ycoord_up(j)   = RotatedCoord(2) + (j-1)*dy*cos(angleRotation);
        ycoord_down(j) = RotatedCoord(2) - j*dy*cos(angleRotation);
        NumMatPoints = NumMatPoints + 2;
    end

    Coord{i,1} = cat(2, xcoord_up, xcoord_down);
    Coord{i,2} = cat(2, ycoord_up, ycoord_down);

    % Store envelope
    x_envelope_upper(i) = xcoord_up(end);
    y_envelope_upper(i) = ycoord_up(end);
    x_envelope_lower(i) = xcoord_down(end);
    y_envelope_lower(i) = ycoord_down(end);

    % Calculate actual height from discretization
    if NumPointsInHeight > 0
        height_actual(i) = (NumPointsInHeight - 1) * dy;
    else
        height_actual(i) = 0;
    end
    height_theoretical(i) = height;

    % Cross-sectional area
    area_cross_section(i) = 2 * height;
end

fprintf('  Generated %d material points in %d sections\n\n', NumMatPoints, BodyNx);

%% ================================================================
%  PART 2: AREA CALCULATIONS
%% ================================================================
fprintf('================================================================\n');
fprintf('  AREA CALCULATIONS\n');
fprintf('================================================================\n\n');

% Method 1: Polygon area from envelope
x_poly = [x_envelope_upper, fliplr(x_envelope_lower)];
y_poly = [y_envelope_upper, fliplr(y_envelope_lower)];
A_polygon = abs(polyarea(x_poly, y_poly));

% Method 2: Trapezoidal integration of cross-sections
x_sections = (0:BodyNx-1) * dx;
A_trapz = trapz(x_sections, area_cross_section);

% Method 3: Analytical approximation
% Head section: integral of sqrt(2*Wh*x - x^2) from 0 to Xh
% Using substitution: integral of semi-circle area formula
% For head (elliptical): A_head ≈ (π/4) * Wh * Xh (approximation)
A_head_analytical = (pi/4) * Wh * Xh;

% Tail section: linear taper
% integral of Wh*(L-x)/(L-Xh) from Xh to L
A_tail_analytical = Wh * (L - Xh) / 2;  % Triangle area

A_analytical_total = A_head_analytical + A_tail_analytical;

% Method 4: Direct summation of all material points (area per point)
area_per_point = dx * dy;
A_material_points = NumMatPoints * area_per_point;

fprintf('1. TOTAL AREA CALCULATIONS:\n');
fprintf('   Method 1 (Polygon):           %.10e\n', A_polygon);
fprintf('   Method 2 (Trapz integration): %.10e\n', A_trapz);
fprintf('   Method 3 (Analytical approx): %.10e\n', A_analytical_total);
fprintf('   Method 4 (Material points):   %.10e\n', A_material_points);
fprintf('\n');
fprintf('   Polygon vs Trapz error:       %.6e (%.4f%%)\n', ...
    abs(A_polygon - A_trapz), 100*abs(A_polygon - A_trapz)/A_trapz);
fprintf('   Polygon vs Analytical error:  %.6e (%.4f%%)\n', ...
    abs(A_polygon - A_analytical_total), 100*abs(A_polygon - A_analytical_total)/A_analytical_total);
fprintf('\n');

%% ================================================================
%  PART 3: VOLUME CALCULATIONS (2D Area × Unit Depth)
%% ================================================================
fprintf('================================================================\n');
fprintf('  VOLUME CALCULATIONS (2D × Unit Depth)\n');
fprintf('================================================================\n\n');

unit_depth = 1.0;  % Assume unit depth in z-direction

V_polygon = A_polygon * unit_depth;
V_trapz = A_trapz * unit_depth;
V_analytical = A_analytical_total * unit_depth;
V_material = A_material_points * unit_depth;

fprintf('Assuming unit depth (z = %.1f):\n', unit_depth);
fprintf('   Volume (Polygon method):      %.10e\n', V_polygon);
fprintf('   Volume (Trapz method):        %.10e\n', V_trapz);
fprintf('   Volume (Analytical):          %.10e\n', V_analytical);
fprintf('   Volume (Material points):     %.10e\n', V_material);
fprintf('\n');

%% ================================================================
%  PART 4: SECTION-WISE ANALYSIS
%% ================================================================
fprintf('================================================================\n');
fprintf('  SECTION-WISE ANALYSIS\n');
fprintf('================================================================\n\n');

% Head section area
A_head_trapz = trapz(x_sections(1:HeadNx), area_cross_section(1:HeadNx));

% Tail section area
if BodyNx > HeadNx
    A_tail_trapz = trapz(x_sections(HeadNx+1:BodyNx), area_cross_section(HeadNx+1:BodyNx));
else
    A_tail_trapz = 0;
end

fprintf('Head Section (0 to %.4f):\n', Xh);
fprintf('   Area (numerical):    %.10e\n', A_head_trapz);
fprintf('   Area (analytical):   %.10e\n', A_head_analytical);
fprintf('   Error:               %.6e (%.4f%%)\n', ...
    abs(A_head_trapz - A_head_analytical), 100*abs(A_head_trapz - A_head_analytical)/A_head_analytical);
fprintf('   Percentage of total: %.2f%%\n', 100*A_head_trapz/A_trapz);
fprintf('\n');

fprintf('Tail Section (%.4f to %.2f):\n', Xh, L);
fprintf('   Area (numerical):    %.10e\n', A_tail_trapz);
fprintf('   Area (analytical):   %.10e\n', A_tail_analytical);
fprintf('   Error:               %.6e (%.4f%%)\n', ...
    abs(A_tail_trapz - A_tail_analytical), 100*abs(A_tail_trapz - A_tail_analytical)/A_tail_analytical);
fprintf('   Percentage of total: %.2f%%\n', 100*A_tail_trapz/A_trapz);
fprintf('\n');

fprintf('Conservation Check:\n');
fprintf('   Head + Tail:         %.10e\n', A_head_trapz + A_tail_trapz);
fprintf('   Total (trapz):       %.10e\n', A_trapz);
fprintf('   Difference:          %.6e\n', abs((A_head_trapz + A_tail_trapz) - A_trapz));
fprintf('\n');

%% ================================================================
%  PART 5: HEIGHT VERIFICATION
%% ================================================================
fprintf('================================================================\n');
fprintf('  CROSS-SECTIONAL HEIGHT VERIFICATION\n');
fprintf('================================================================\n\n');

% Calculate errors in height discretization
height_error = height_actual - height_theoretical;
height_error_rel = abs(height_error) ./ (height_theoretical + eps);

max_height_error = max(abs(height_error));
mean_height_error = mean(abs(height_error));
max_height_error_rel = max(height_error_rel(height_theoretical > 0));

fprintf('Height Discretization Analysis:\n');
fprintf('   Max height (theoretical): %.6e at x = %.4f\n', ...
    max(height_theoretical), x_sections(height_theoretical == max(height_theoretical)));
fprintf('   Max height (actual):      %.6e\n', max(height_actual));
fprintf('   Max absolute error:       %.6e\n', max_height_error);
fprintf('   Mean absolute error:      %.6e\n', mean_height_error);
fprintf('   Max relative error:       %.4f%%\n', 100*max_height_error_rel);
fprintf('\n');

%% ================================================================
%  PART 6: PERIMETER CALCULATION
%% ================================================================
fprintf('================================================================\n');
fprintf('  PERIMETER CALCULATION\n');
fprintf('================================================================\n\n');

% Calculate perimeter from envelope
dx_upper = diff(x_envelope_upper);
dy_upper = diff(y_envelope_upper);
ds_upper = sqrt(dx_upper.^2 + dy_upper.^2);
L_upper = sum(ds_upper);

dx_lower = diff(x_envelope_lower);
dy_lower = diff(y_envelope_lower);
ds_lower = sqrt(dx_lower.^2 + dy_lower.^2);
L_lower = sum(ds_lower);

% Closing segments (nose and tail)
dx_nose = x_envelope_lower(1) - x_envelope_upper(1);
dy_nose = y_envelope_lower(1) - y_envelope_upper(1);
L_nose = sqrt(dx_nose^2 + dy_nose^2);

dx_tail = x_envelope_upper(end) - x_envelope_lower(end);
dy_tail = y_envelope_upper(end) - y_envelope_lower(end);
L_tail = sqrt(dx_tail^2 + dy_tail^2);

Perimeter = L_upper + L_lower + L_nose + L_tail;

fprintf('Perimeter Components:\n');
fprintf('   Upper surface:  %.6e\n', L_upper);
fprintf('   Lower surface:  %.6e\n', L_lower);
fprintf('   Nose closure:   %.6e\n', L_nose);
fprintf('   Tail closure:   %.6e\n', L_tail);
fprintf('   Total:          %.6e\n', Perimeter);
fprintf('\n');

% Compactness metric (isoperimetric ratio)
% For a circle: 4πA/P² = 1
% For other shapes: < 1
compactness = 4 * pi * A_polygon / (Perimeter^2);
fprintf('Shape Metrics:\n');
fprintf('   Compactness (4πA/P²): %.6f\n', compactness);
fprintf('   Aspect ratio (L/Wh):  %.2f\n', L/Wh);
fprintf('\n');

%% ================================================================
%  PART 7: VISUALIZATION
%% ================================================================
fprintf('Generating visualization...\n');

fig = figure('Visible', 'on', 'Position', [100 100 1800 1200]);

% Subplot 1: Full geometry with envelope
subplot(2,3,1);
hold on;
% Plot all material points
for i = 1:size(Coord,1)
    plot(Coord{i,1}(:), Coord{i,2}(:), '.', 'MarkerSize', 4, 'Color', [0.7 0.7 0.7]);
end
% Plot envelope
plot(x_envelope_upper, y_envelope_upper, '-', 'LineWidth', 2.5, 'Color', [0.85 0.15 0.15], 'DisplayName', 'Upper envelope');
plot(x_envelope_lower, y_envelope_lower, '-', 'LineWidth', 2.5, 'Color', [0.1 0.4 0.8], 'DisplayName', 'Lower envelope');
% Mark head-tail boundary
xline(Xh - 0.5, '--k', 'LineWidth', 1.5, 'DisplayName', 'Head/Tail boundary');
axis equal; grid on; box on;
xlabel('x', 'FontSize', 12); ylabel('y', 'FontSize', 12);
title('Eel Geometry with Envelope', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
set(gca, 'FontSize', 11);

% Subplot 2: Cross-sectional area distribution
subplot(2,3,2);
area(x_sections, area_cross_section, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', [0.1 0.3 0.6], 'LineWidth', 1.5);
hold on;
xline(Xh, '--r', 'LineWidth', 2, 'DisplayName', 'Head/Tail boundary');
grid on; box on;
xlabel('x [m]', 'FontSize', 12);
ylabel('Cross-sectional Area [m]', 'FontSize', 12);
title('Cross-Sectional Area Distribution', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);

% Subplot 3: Height verification
subplot(2,3,3);
plot(x_sections, height_theoretical, '-', 'LineWidth', 2.5, 'Color', [0.85 0.15 0.15], 'DisplayName', 'Theoretical');
hold on;
plot(x_sections, height_actual, '--', 'LineWidth', 2.0, 'Color', [0.1 0.4 0.8], 'DisplayName', 'Discretized');
xline(Xh, ':k', 'LineWidth', 1.5);
grid on; box on;
xlabel('x [m]', 'FontSize', 12);
ylabel('Half-Height [m]', 'FontSize', 12);
title('Cross-Sectional Height', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);

% Subplot 4: Area comparison bar chart
subplot(2,3,4);
bar_data = [A_polygon, A_trapz, A_analytical_total];
bar_labels = {'Polygon', 'Trapz', 'Analytical'};
bar(bar_data, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
set(gca, 'XTickLabel', bar_labels);
ylabel('Area [m^2]', 'FontSize', 12);
title('Total Area Comparison', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
for i = 1:length(bar_data)
    text(i, bar_data(i)*1.05, sprintf('%.6f', bar_data(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end
set(gca, 'FontSize', 11);

% Subplot 5: Section-wise area comparison
subplot(2,3,5);
section_data = [A_head_trapz, A_tail_trapz; A_head_analytical, A_tail_analytical];
bar(section_data', 'grouped');
set(gca, 'XTickLabel', {'Head', 'Tail'});
ylabel('Area [m^2]', 'FontSize', 12);
title('Section-Wise Area', 'FontSize', 14, 'FontWeight', 'bold');
legend('Numerical', 'Analytical', 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);

% Subplot 6: Height error distribution
subplot(2,3,6);
plot(x_sections, abs(height_error), '-o', 'LineWidth', 2.0, 'MarkerSize', 4, 'Color', [0.8 0.3 0.3]);
hold on;
xline(Xh, '--k', 'LineWidth', 1.5);
grid on; box on;
xlabel('x [m]', 'FontSize', 12);
ylabel('|Height Error| [m]', 'FontSize', 12);
title('Height Discretization Error', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);

sgtitle('Eel2D Volume and Area Verification', 'FontSize', 18, 'FontWeight', 'bold');

% Save figure
if ~exist('output', 'dir'), mkdir('output'); end
if ~exist('output/images', 'dir'), mkdir('output/images'); end
save_file = fullfile('output', 'images', 'volume_area_verification.png');
saveas(fig, save_file);
fprintf('Saved: %s\n\n', save_file);

%% ================================================================
%  PART 8: SAVE DATA TO FILE
%% ================================================================
fprintf('Saving verification data...\n');

if ~exist('output/data', 'dir'), mkdir('output/data'); end
data_file = fullfile('output', 'data', 'volume_area_data.txt');
fid = fopen(data_file, 'wt');

fprintf(fid, '================================================================\n');
fprintf(fid, '  EEL2D VOLUME AND AREA VERIFICATION DATA\n');
fprintf(fid, '================================================================\n\n');

fprintf(fid, 'CONFIGURATION:\n');
fprintf(fid, '  Body length (L):     %.6e m\n', L);
fprintf(fid, '  Head width (Wh):     %.6e m\n', Wh);
fprintf(fid, '  Head length (Xh):    %.6e m\n', Xh);
fprintf(fid, '  Grid spacing (dx):   %.6e m\n', dx);
fprintf(fid, '  Grid spacing (dy):   %.6e m\n', dy);
fprintf(fid, '  Total sections:      %d\n', BodyNx);
fprintf(fid, '  Material points:     %d\n\n', NumMatPoints);

fprintf(fid, 'AREA CALCULATIONS:\n');
fprintf(fid, '  Polygon method:      %.10e m^2\n', A_polygon);
fprintf(fid, '  Trapz method:        %.10e m^2\n', A_trapz);
fprintf(fid, '  Analytical approx:   %.10e m^2\n', A_analytical_total);
fprintf(fid, '  Material points:     %.10e m^2\n\n', A_material_points);

fprintf(fid, 'VOLUME CALCULATIONS (unit depth):\n');
fprintf(fid, '  Polygon method:      %.10e m^3\n', V_polygon);
fprintf(fid, '  Trapz method:        %.10e m^3\n', V_trapz);
fprintf(fid, '  Analytical approx:   %.10e m^3\n\n', V_analytical);

fprintf(fid, 'SECTION-WISE AREAS:\n');
fprintf(fid, '  Head (numerical):    %.10e m^2\n', A_head_trapz);
fprintf(fid, '  Head (analytical):   %.10e m^2\n', A_head_analytical);
fprintf(fid, '  Tail (numerical):    %.10e m^2\n', A_tail_trapz);
fprintf(fid, '  Tail (analytical):   %.10e m^2\n\n', A_tail_analytical);

fprintf(fid, 'PERIMETER:\n');
fprintf(fid, '  Total perimeter:     %.10e m\n', Perimeter);
fprintf(fid, '  Compactness:         %.6f\n', compactness);
fprintf(fid, '  Aspect ratio:        %.2f\n\n', L/Wh);

fprintf(fid, 'HEIGHT VERIFICATION:\n');
fprintf(fid, '  Max height:          %.6e m\n', max(height_theoretical));
fprintf(fid, '  Max abs error:       %.6e m\n', max_height_error);
fprintf(fid, '  Mean abs error:      %.6e m\n', mean_height_error);
fprintf(fid, '  Max rel error:       %.4f%%\n\n', 100*max_height_error_rel);

fclose(fid);
fprintf('Saved: %s\n\n', data_file);

%% ================================================================
%  SUMMARY
%% ================================================================
fprintf('================================================================\n');
fprintf('  VERIFICATION SUMMARY\n');
fprintf('================================================================\n\n');

fprintf('✓ Total Area:        %.6e m^2\n', A_polygon);
fprintf('✓ Total Volume:      %.6e m^3 (unit depth)\n', V_polygon);
fprintf('✓ Head Area:         %.6e m^2 (%.1f%% of total)\n', A_head_trapz, 100*A_head_trapz/A_trapz);
fprintf('✓ Tail Area:         %.6e m^2 (%.1f%% of total)\n', A_tail_trapz, 100*A_tail_trapz/A_trapz);
fprintf('✓ Perimeter:         %.6e m\n', Perimeter);
fprintf('✓ Material Points:   %d\n', NumMatPoints);
fprintf('✓ Aspect Ratio:      %.1f:1\n', L/Wh);
fprintf('\n');

fprintf('Accuracy Metrics:\n');
fprintf('  Area methods agree within:    %.4f%%\n', 100*abs(A_polygon - A_trapz)/A_trapz);
fprintf('  Analytical error:             %.4f%%\n', 100*abs(A_polygon - A_analytical_total)/A_analytical_total);
fprintf('  Height discretization error:  %.4f%%\n', 100*max_height_error_rel);
fprintf('  Conservation check:           %.6e (< 1e-10 = excellent)\n', abs((A_head_trapz + A_tail_trapz) - A_trapz));
fprintf('\n');

fprintf('================================================================\n');
fprintf('  TEST COMPLETE\n');
fprintf('================================================================\n\n');
