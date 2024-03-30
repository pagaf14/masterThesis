clear;
clc;
close all;

% OSSERVAZIONI
%
% La sensitività superficiale media aumenta (?)
% DSN_001 --> DNS_020, 10 iterazioni, da DSN_=010 solo cartella DEFORM
% Legame iterazioni - DSN?

%% Geometry

surface_flow_1DIR=readtable("../DESIGNS/DSN_001/DIRECT/surface_flow.csv");
x_points_1DIR=surface_flow_1DIR.x;
y_points_1DIR=surface_flow_1DIR.y;

surface_flow_1AD=readtable("../DESIGNS/DSN_001/ADJOINT_DRAG/surface_flow.csv");
x_points_1ADf=surface_flow_1AD.x;
y_points_1ADf=surface_flow_1AD.y;

surface_adjoint_1AD=readtable("../DESIGNS/DSN_001/ADJOINT_DRAG/surface_adjoint.csv");
x_points_1ADadj=surface_adjoint_1AD.x;
y_points_1ADadj=surface_adjoint_1AD.y;

figure;
hold on;
plot(x_points_1DIR, y_points_1DIR);
plot(x_points_1ADf, y_points_1ADf);
plot(x_points_1ADadj, y_points_1ADadj);
legend('DIRECT', 'SURFACE_FLOW', 'SURFACE_ADJOINT');

% ==> i punti in questi 3 file sono coincidenti!!

%% Pressure coefficient

surface_flow_1=readtable("../DESIGNS/DSN_001/DIRECT/surface_flow.csv");
cp_1=surface_flow_1.Pressure_Coefficient;
x_points_1=surface_flow_1.x;
y_points_1=surface_flow_1.y;

% surface_flow_2=readtable("../DESIGNS/DSN_002/DIRECT/surface_flow.csv");
% cp_2=surface_flow_2.Pressure_Coefficient;
% x_points_2=surface_flow_2.x;
% y_points_2=surface_flow_2.y;

surface_flow_20=readtable("../DESIGNS/DSN_020/DIRECT/surface_flow.csv");
cp_20=surface_flow_20.Pressure_Coefficient;
x_points_20=surface_flow_20.x;
y_points_20=surface_flow_20.y;

figure;
hold on;
plot1=plot(x_points_1, -3*y_points_1+1.5, 'Color', '#88ACE0', 'LineWidth', 0.5);
plot(x_points_1, cp_1,'Color', '#88ACE0');

% plot2=plot(x_points_2, -3*y_points_2+1.5, 'Color', '#199FD6', 'LineWidth', 0.5);
% plot(x_points_2, cp_2, 'Color', '#199FD6');

plot2=plot(x_points_20, -3*y_points_20+1.5, 'Color', '#0047AB', 'LineWidth', 0.5);
plot(x_points_20, cp_20, 'Color', '#0047AB');

legend([plot1, plot2], 'NACA $0012$', 'OPTIMIZED');
xlabel('$x/c$');
ylabel('$C_P$');
set(gca, 'YDir', 'reverse');

exportStandardizedFigure(gcf, 'NACA0012_CP', 1, 'legendLocation', 'NorthEast','legendOrientation', 'vertical','addMarkers', false, 'changeColors', false);


%% Drag history

history_project=readtable("../history_project.csv");
history_eval=history_project.EVALUATION;
history_drag=history_project.DRAG;

figure;
plot(history_eval, history_drag, 'Color', '#C4151C', 'Marker', '*');
xlabel('EVALUATION');
ylabel('DRAG');

exportStandardizedFigure(gcf, 'NACA0012_HISTORY_DRAG', 1, 'legendLocation', 'NorthEast','legendOrientation', 'vertical', 'addMarkers', false, 'changeColors', false);

%% Mean surface sensitivity history

num_folders = 9;  % Numero totale di cartelle
surface_sensitivity = zeros(1, num_folders);  % Inizializza il vettore dei risultati

for i = 1:num_folders
    folder_name = sprintf('../DESIGNS/DSN_%03d/', i);  % Genera il nome della cartella dinamicamente
    file_path = fullfile(folder_name, 'ADJOINT_DRAG', 'surface_flow.csv');  % Genera il percorso del file
    surface_flow = readtable(file_path);  % Legge il file
    surface_sensitivity(i) = mean(surface_flow.Surface_Sensitivity);  % Calcola la media della colonna Surface_Sensitivity
end

figure;
plot([1:num_folders], surface_sensitivity, 'Color', '#FF8000', 'Marker', '*');
xlabel('Iterations');
ylabel('Mean surface sensitivity')

exportStandardizedFigure(gcf, 'NACA0012_HISTORY_SURF_SENS', 1, 'legendLocation', 'NorthEast','legendOrientation', 'vertical', 'addMarkers', false, 'changeColors', false);
