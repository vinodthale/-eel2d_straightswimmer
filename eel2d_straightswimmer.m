%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2014 - 2019 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

% Filename eel2D_straightswimmer.m
% Created By Amneet Bhalla

clear all;
clc;

% Create output directories if they don't exist
if ~exist('output', 'dir'), mkdir('output'); end
if ~exist('output/vertices', 'dir'), mkdir('output/vertices'); end
if ~exist('output/images', 'dir'), mkdir('output/images'); end
if ~exist('output/pdfs', 'dir'), mkdir('output/pdfs'); end
if ~exist('output/data', 'dir'), mkdir('output/data'); end

%Background Mesh Properties
Ny = 16*4*4*4; Ly = 4;
Nx = 32*4*4*4; Lx = 8;
dy = Ly/Ny; dx = Lx/Nx;

L  = 1;                 %length of the body
Wh = 0.04*L;            %width of the head.
Xh = 0.04;   

% No of points on backbone.
BodyNx = ceil(L/dx);
HeadNx = ceil(0.04/dx);
Coord{BodyNx,2} = [];
NumMatPoints    =  0;

angleRotation  = 0;
RotationMatrix = [ cos(angleRotation) -sin(angleRotation);
                   sin(angleRotation)  cos(angleRotation)];

%store from nose to head of the body.
for i = 1:HeadNx
     
    x                 = (i-1)*dx;
    y                 = 0.125* ((x + 0.03125)/1.03125)*sin(2*pi*x);
    height            = sqrt(2*Wh*x - x^2);
    NumPointsInHeight = ceil(height/dy);
    
    ycoord_up = []; ycoord_down = []; xcoord_up = []; xcoord_down = [];
    for j = 1:NumPointsInHeight
       
        xshifted        = x - 0.5;
        yshifted        = y;
        RotatedCoord    = RotationMatrix*[xshifted;yshifted];
        
        
        xcoord_up(j)   = RotatedCoord(1) -(j-1)*dy*sin(angleRotation);
        xcoord_down(j) = RotatedCoord(1) + j*dy*sin(angleRotation);
        ycoord_up(j)   = RotatedCoord(2) + (j-1)*dy*cos(angleRotation);
        ycoord_down(j) = RotatedCoord(2) -j*dy*cos(angleRotation);
        NumMatPoints   = NumMatPoints+2;
    end
    
    Coord{i,1} = cat(2,xcoord_up,xcoord_down);
    Coord{i,2} = cat(2,ycoord_up,ycoord_down);
    
end


%store from head to tail of the body.
for i = HeadNx+1 : BodyNx
    
    x     = (i-1)*dx;
    y     = 0.125* ((x + 0.03125)/1.03125)*sin(2*pi*x);
    height            = Wh*(L - x)/(L - Xh);
    NumPointsInHeight = ceil(height/dy);
    
    ycoord_up = []; ycoord_down = []; xcoord_up = []; xcoord_down = [];
    for j = 1:NumPointsInHeight
        
        xshifted        = x - 0.5;
        yshifted        = y;
        RotatedCoord    = RotationMatrix*[xshifted;yshifted];
        
        xcoord_up(j)   = RotatedCoord(1) -(j-1)*dy*sin(angleRotation);
        xcoord_down(j) = RotatedCoord(1) + j*dy*sin(angleRotation);
        ycoord_up(j)   = RotatedCoord(2) + (j-1)*dy*cos(angleRotation);
        ycoord_down(j) = RotatedCoord(2) -j*dy*cos(angleRotation);
        NumMatPoints = NumMatPoints+2;
    end
    
    Coord{i,1} = cat(2,xcoord_up,xcoord_down);
    Coord{i,2} = cat(2,ycoord_up,ycoord_down);
    
end
    

%% Calculate geometric properties
fprintf('Calculating geometric properties...\n');

% Extract envelope (upper and lower surfaces)
x_envelope_upper = zeros(1, size(Coord,1));
y_envelope_upper = zeros(1, size(Coord,1));
x_envelope_lower = zeros(1, size(Coord,1));
y_envelope_lower = zeros(1, size(Coord,1));

for i = 1:size(Coord,1)
    if ~isempty(Coord{i,1}) && ~isempty(Coord{i,2})
        x_coords = Coord{i,1};
        y_coords = Coord{i,2};

        % Find uppermost and lowermost points
        [y_max, idx_max] = max(y_coords);
        [y_min, idx_min] = min(y_coords);

        x_envelope_upper(i) = x_coords(idx_max);
        y_envelope_upper(i) = y_max;
        x_envelope_lower(i) = x_coords(idx_min);
        y_envelope_lower(i) = y_min;
    end
end

% Calculate total area using polygon method
x_poly = [x_envelope_upper, fliplr(x_envelope_lower)];
y_poly = [y_envelope_upper, fliplr(y_envelope_lower)];
total_area = abs(polyarea(x_poly, y_poly));

% Calculate volume (assuming unit depth in z-direction)
unit_depth = 1.0;
total_volume = total_area * unit_depth;

% Calculate perimeter
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

perimeter = L_upper + L_lower + L_nose + L_tail;

% Compactness metric (4πA/P² - equals 1 for circle, < 1 for other shapes)
compactness = 4 * pi * total_area / (perimeter^2);

fprintf('\n');
fprintf('================================================================\n');
fprintf('  GEOMETRIC PROPERTIES\n');
fprintf('================================================================\n');
fprintf('Body dimensions:\n');
fprintf('  Length: %.4f m\n', L);
fprintf('  Max width: %.4f m\n', Wh);
fprintf('\n');
fprintf('Calculated properties:\n');
fprintf('  Total area:    %.6e m^2\n', total_area);
fprintf('  Total volume:  %.6e m^3 (unit depth)\n', total_volume);
fprintf('  Perimeter:     %.6e m\n', perimeter);
fprintf('  Compactness:   %.4f (1.0 = circle)\n', compactness);
fprintf('  Aspect ratio:  %.1f:1\n', L/Wh);
fprintf('================================================================\n\n');

% plot the fish
for i = 1:size(Coord,1)

    plot(Coord{i,1}(:),Coord{i,2}(:),'.');
    hold on;
end

% Plot envelope
plot(x_envelope_upper, y_envelope_upper, '-r', 'LineWidth', 2, 'DisplayName', 'Upper envelope');
plot(x_envelope_lower, y_envelope_lower, '-b', 'LineWidth', 2, 'DisplayName', 'Lower envelope');
legend('show');


% write the coordinates into txt file.
vertex_file = fullfile('output', 'vertices', 'eel2d.vertex');
fid = fopen(vertex_file,'wt');
fprintf(fid,'%d\n',NumMatPoints);

for i = 1:size(Coord,1)

    for j = 1: length( Coord{i,1}  )

         fprintf(fid,'%f\t%f\n', Coord{i,1}(j), Coord{i,2}(j) );
    end

end
fclose(fid);

% Save the plot
fig = gcf;
set(fig, 'Position', [100 100 1200 800]);
xlabel('x', 'FontSize', 12);
ylabel('y', 'FontSize', 12);
title('Eel 2D Straight Swimmer', 'FontSize', 14, 'FontWeight', 'bold');
axis equal;
grid on;
plot_file = fullfile('output', 'images', 'eel2d_straightswimmer.png');
saveas(fig, plot_file);
fprintf('Saved vertex file: %s\n', vertex_file);
fprintf('Saved plot: %s\n', plot_file);

% Save geometric properties to data file
data_file = fullfile('output', 'data', 'eel_geometry_properties.txt');
fid_data = fopen(data_file, 'wt');
fprintf(fid_data, '================================================================\n');
fprintf(fid_data, '  EEL2D GEOMETRIC PROPERTIES\n');
fprintf(fid_data, '================================================================\n\n');
fprintf(fid_data, 'Body dimensions:\n');
fprintf(fid_data, '  Length (L):           %.6e m\n', L);
fprintf(fid_data, '  Head width (Wh):      %.6e m\n', Wh);
fprintf(fid_data, '  Head length (Xh):     %.6e m\n', Xh);
fprintf(fid_data, '  Grid spacing (dx):    %.6e m\n', dx);
fprintf(fid_data, '  Grid spacing (dy):    %.6e m\n\n', dy);
fprintf(fid_data, 'Material points:\n');
fprintf(fid_data, '  Total points:         %d\n', NumMatPoints);
fprintf(fid_data, '  Sections:             %d\n\n', BodyNx);
fprintf(fid_data, 'Calculated properties:\n');
fprintf(fid_data, '  Total area:           %.10e m^2\n', total_area);
fprintf(fid_data, '  Total volume:         %.10e m^3 (unit depth)\n', total_volume);
fprintf(fid_data, '  Perimeter:            %.10e m\n', perimeter);
fprintf(fid_data, '  Upper surface length: %.10e m\n', L_upper);
fprintf(fid_data, '  Lower surface length: %.10e m\n', L_lower);
fprintf(fid_data, '  Compactness (4πA/P²): %.6f\n', compactness);
fprintf(fid_data, '  Aspect ratio (L/Wh):  %.2f:1\n\n', L/Wh);
fprintf(fid_data, '================================================================\n');
fclose(fid_data);
fprintf('Saved geometry data: %s\n', data_file);
