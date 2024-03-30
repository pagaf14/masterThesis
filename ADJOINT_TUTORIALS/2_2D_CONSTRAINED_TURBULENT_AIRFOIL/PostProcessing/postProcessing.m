clear;
clc;
close all;

%% Geometry

surface_flow_1DIR=readtable("../DESIGNS/DSN_001/DIRECT/surface_flow.csv");
x_points_1DIR=surface_flow_1DIR.x;
y_points_1DIR=surface_flow_1DIR.y;

surface_flow_2DIR=readtable("../DESIGNS/DSN_002/DIRECT/surface_flow.csv");
x_points_2DIR=surface_flow_2DIR.x;
y_points_2DIR=surface_flow_2DIR.y;

surface_flow_2DEF=readtable("../DESIGNS/DSN_002/DEFORM/surface_deformed.csv");
x_points_2DEF=surface_flow_2DEF.x;
y_points_2DEF=surface_flow_2DEF.y;

figure;
hold on;
plot(x_points_1DIR, y_points_1DIR);
plot(x_points_2DIR, y_points_2DIR);
plot(x_points_2DEF, y_points_2DEF);
xlim([0 1]);
ylim([-0.5 0.5]);
legend('DIRECT_1', 'DIRECT_2', 'DEFORM_2');

% Tra DSN_001 e DSN_002 la geometria cambia, ma in un design tra DIRECT e
% DEFORM non cambia

%% Pressure coefficient

surface_flow_1=readtable("../DESIGNS/DSN_001/DIRECT/surface_flow.csv");
cp_1=surface_flow_1.Pressure_Coefficient;
x_points_1=surface_flow_1.x;
y_points_1=surface_flow_1.y;

surface_flow_14=readtable("../DESIGNS/DSN_014/DIRECT/surface_flow.csv");
cp_14=surface_flow_14.Pressure_Coefficient;
x_points_14=surface_flow_14.x;
y_points_14=surface_flow_14.y;

figure;
hold on;
plot1=plot(x_points_1, -3*y_points_1+1.5, '-k', 'LineWidth', 0.05);
plot(x_points_1, cp_1,'-k');

plot2=plot(x_points_14, -3*y_points_14+1.5, 'Color', '#C4151C', 'LineWidth', 0.05);
plot(x_points_14, cp_14, 'Color', '#C4151C');

legend([plot1, plot2], 'RAE $2822$', 'OPTIMIZED AIRFOIL')
xlabel('$x/c$');
ylabel('$C_P$');
set(gca, 'YDir', 'reverse');

exportStandardizedFigure(gcf, 'RAE2822_CP', 1, 'legendLocation', 'NorthEast','legendOrientation', 'vertical','addMarkers', false, 'changeColors', false, 'changeLineStyle', false);

% Il grafico non si raccorda al bordo d'uscita

%% Drag and lift history

history_project=readtable("../history_project.csv");
history_eval=history_project.EVALUATION;
history_drag=history_project.DRAG;
history_lift=history_project.LIFT;

figure;
hold on;

yyaxis left; % Utilizza l'asse y sinistro
plot(history_eval, history_drag, 'Color', '#C4151C', 'Marker', '*', 'LineWidth', 1.5);
ylabel('C_D');
ylim([min(history_drag), max(history_drag)]);
set(gca, 'YColor', 'k');

yyaxis right; % Utilizza l'asse y destro
plot(history_eval, history_lift, 'Color', '#0047AB', 'Marker', '*', 'LineWidth', 1.5);
ylabel('$C_L$');
ylim([0 0.8]);
set(gca, 'YColor', 'k');

xlabel('EVALUATION');
legend('$C_D$', '$C_L$')

exportStandardizedFigure(gcf, 'RAE2822_HISTORY_CD_CL', 1, 'legendLocation', 'SouthWest','legendOrientation', 'vertical', 'addMarkers', false, 'changeColors', false, 'changeLineStyle', false);

%% Pitching moment and airfoil thickness history

history_moment_z=history_project.MOMENT_Z;
history_airfoil_thickness=history_project.AIRFOIL_THICKNESS;

figure;
hold on;

yyaxis left; % Utilizza l'asse y sinistro
plot(history_eval, history_moment_z, 'Color', '#FF8000', 'Marker', '*', 'LineWidth', 1.5);
ylabel('MOMENT\_Z');
ylim([0.07 0.095]);
set(gca, 'YColor', 'k');

yyaxis right; % Utilizza l'asse y destro
plot(history_eval, history_airfoil_thickness, 'Color', '#FDAF00', 'Marker', '*', 'LineWidth', 1.5);
ylabel('AIRFOIL_THICKNESS');
ylim([0.12 0.14]);
set(gca, 'YColor', 'k');

xlabel('EVALUATION');
legend('MOMENT_Z', 'AIRFOIL_THICKNESS')

exportStandardizedFigure(gcf, 'RAE2822_HISTORY_MZ_AT', 1, 'legendLocation', 'NorthEast','legendOrientation', 'vertical', 'addMarkers', false, 'changeColors', false, 'changeLineStyle', false);
