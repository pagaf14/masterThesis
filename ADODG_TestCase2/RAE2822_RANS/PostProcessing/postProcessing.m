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

surface_flow_119DIR=readtable("../DESIGNS/DSN_119/DIRECT/surface_flow.csv");
x_points_119DIR=surface_flow_119DIR.x;
y_points_119DIR=surface_flow_119DIR.y;

paper_geometry = readtable("OPTIMIZED_GEOMETRY_paper.csv");
x_points_paper = paper_geometry.Var1;
y_points_paper = paper_geometry.Var2;

figure;
hold on;
plot([x_points_1DIR; x_points_1DIR(1)], [y_points_1DIR; y_points_1DIR(1)], '-k', 'LineWidth', 0.05);
% plot(x_points_2DIR, y_points_2DIR);
% plot(x_points_2DEF, y_points_2DEF);
plot([x_points_119DIR; x_points_119DIR(1)], [y_points_119DIR; y_points_119DIR(1)], 'Color', '#C4151C', 'LineWidth', 0.05);
plot([x_points_paper; x_points_paper(1)], [y_points_paper; y_points_paper(1)], '--', 'Color', '#C4151C', 'LineWidth', 0.05);
xlim([0 1]);
ylim([-0.1 0.08]);
xlabel("x/c");
ylabel("z");
legend('BASELINE GEOMETRY', 'OPTIMIZED GEOMETRY (SU2)', 'OPTIMIZED GEOMETRY (AMRIT 2017)');

exportStandardizedFigure(gcf, 'RAE2822_GEOMETRY', 1, 'legendLocation', 'SouthEast','legendOrientation', 'vertical','addMarkers', false, 'changeColors', false, 'changeLineStyle', false, 'overwriteFigure', true);

% Tra DSN_001 e DSN_002 la geometria cambia, ma in un design tra DIRECT e
% DEFORM non cambia

%% Pressure coefficient

surface_flow_1=readtable("../DESIGNS/DSN_001/DIRECT/surface_flow.csv");
cp_1=surface_flow_1.Pressure_Coefficient;
x_points_1=surface_flow_1.x;
y_points_1=surface_flow_1.y;

surface_flow_119=readtable("../DESIGNS/DSN_119/DIRECT/surface_flow.csv");
cp_119=surface_flow_119.Pressure_Coefficient;
x_points_119=surface_flow_119.x;
y_points_119=surface_flow_119.y;

paper_geometry = readtable("OPTIMIZED_GEOMETRY_paper.csv");
x_points_paper = paper_geometry.Var1;
y_points_paper = paper_geometry.Var2;

paper_cp = readtable("CP_paper.csv");
cp_paper = paper_cp.Var2;
x_points_cp_paper = paper_cp.Var1;

figure;
hold on;
plot1 = plot([x_points_1; x_points_1(1)], [-5*y_points_1+1.5; -5*y_points_1(1)+1.5], '-k', 'LineWidth', 0.05);
plot([x_points_1; x_points_1(1)], [cp_1; cp_1(1)],'-k');

plot2 = plot([x_points_119; x_points_119(1)], [-5*y_points_119+1.5; -5*y_points_119(1)+1.5], 'Color', '#C4151C', 'LineWidth', 0.05);
plot([x_points_119; x_points_119(1)], [cp_119; cp_119(1)], 'Color', '#C4151C');

plot3 = plot([x_points_cp_paper; x_points_cp_paper(1)], [-cp_paper; -cp_paper(1)], '--', 'Color', '#C4151C', 'LineWidth', 0.05);
plot([x_points_paper; x_points_paper(1)], [-5*y_points_paper+1.5; -5*y_points_paper(1)+1.5], '--', 'Color', '#C4151C', 'LineWidth', 0.05);

legend([plot1, plot2, plot3], 'RAE $2822$', 'OPTIMIZED AIRFOIL (SU2)', 'OPTIMIZED AIRFOIL (AMRIT 2017)', 'FontSize', 5);
xlim([0 1.1]);
ylim([-1.5 2.5])
xlabel('$x/c$');
ylabel('$C_P$');
set(gca, 'YDir', 'reverse');

exportStandardizedFigure(gcf, 'RAE2822_CP', 1, 'legendLocation', 'SouthEast','legendOrientation', 'vertical','addMarkers', false, 'changeColors', false, 'changeLineStyle', false, 'overwriteFigure', true);

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
ylim([min(history_drag)-0.0005, max(history_drag)+0.001]);
set(gca, 'YColor', 'k');

yyaxis right; % Utilizza l'asse y destro
plot(history_eval, history_lift, 'Color', '#0047AB', 'Marker', '*', 'LineWidth', 1.5);
ylabel('$C_L$');
ylim([0 0.85]);
set(gca, 'YColor', 'k');

xlabel('EVALUATION');
legend('$C_D$', '$C_L$')

exportStandardizedFigure(gcf, 'RAE2822_HISTORY_CD_CL', 1, 'legendLocation', 'best','legendOrientation', 'vertical', 'addMarkers', false, 'changeColors', false, 'changeLineStyle', false, 'overwriteFigure', true);

%% Pitching moment and airfoil area history

history_moment_z=history_project.MOMENT_Z;
history_airfoil_area=history_project.AIRFOIL_AREA;

figure;
hold on;

yyaxis left; % Utilizza l'asse y sinistro
plot(history_eval, -history_moment_z, 'Color', '#FF8000', 'Marker', '*', 'LineWidth', 1.5);
ylabel('C_M');
%ylim([0.07 0.095]);
set(gca, 'YColor', 'k');
hold on;
plot(history_eval, -0.092*ones(size(history_eval)),'--', 'Color', '#FF8000', 'LineWidth', 1.5);
xlim([0 130]);

yyaxis right; % Utilizza l'asse y destro
plot(history_eval, history_airfoil_area, 'Color', '#6CC407', 'Marker', '*', 'LineWidth', 1.5);
ylabel('AREA');
set(gca, 'YColor', 'k');
hold on;
plot(history_eval, 0.0778446*ones(size(history_eval)), '--', 'Color', '#6CC407', 'LineWidth', 1.5);

xlabel('EVALUATION');
legend('C_M', 'C_M = -0.092', 'AIRFOIL AREA', 'BASELINE AIRFOIL AREA')

exportStandardizedFigure(gcf, 'RAE2822_HISTORY_MZ_AA', 1, 'legendLocation', 'NorthEast','legendOrientation', 'vertical', 'addMarkers', false, 'changeColors', false, 'changeLineStyle', false, 'overwriteFigure', true);

%% Surface sensitivities
% ITERATION #1
ss_DRAG_1 = readtable('../DESIGNS/DSN_001/ADJOINT_DRAG/surface_sens.csv');
ss_DRAG_1 = ss_DRAG_1.Surface_Sensitivity;

ss_MOMENT_1 = readtable('../DESIGNS/DSN_001/ADJOINT_MOMENT_Z/surface_adjoint.csv');
ss_MOMENT_1 = ss_MOMENT_1.Surface_Sensitivity;

surface_flow_1=readtable("../DESIGNS/DSN_001/DIRECT/surface_flow.csv");
x_points_1=surface_flow_1.x;

figure;
hold on;
plot(x_points_1, ss_DRAG_1, 'Color', '#B100FF', 'Marker', '*');
plot(x_points_1, ss_MOMENT_1, 'Color', '#88ACE0', 'Marker', '*');
legend('DRAG', 'MOMENT_Z');
xlabel('$x/c$');
ylabel('Surface sensitivities');
ylim([-3 2]);
title('Surface sensitivities ITERATION #1')

file_name = "RAE2822_SURFACE_SENSITIVITIES_ITER1.pdf"
if exist(file_name) == 2
    delete(file_name)
end

exportStandardizedFigure(gcf, 'RAE2822_SURFACE_SENSITIVITIES_ITER1  ', 1, 'legendLocation', 'SouthWest','legendOrientation', 'vertical','addMarkers', false, 'changeColors', false);


% ITERATION #119

ss_DRAG_119 = readtable('../DESIGNS/DSN_119/ADJOINT_DRAG/surface_sens.csv');
ss_DRAG_119 = ss_DRAG_119.Surface_Sensitivity;

ss_MOMENT_119 = readtable('../DESIGNS/DSN_119/ADJOINT_MOMENT_Z/surface_adjoint.csv');
ss_MOMENT_119 = ss_MOMENT_119.Surface_Sensitivity;

surface_flow_119=readtable("../DESIGNS/DSN_119/DIRECT/surface_flow.csv");
x_points_119=surface_flow_119.x;

figure;
hold on;
plot(x_points_119, ss_DRAG_119, 'Color', '#B100FF', 'Marker', '*');
plot(x_points_119, ss_MOMENT_119, 'Color', '#88ACE0', 'Marker', '*');
legend('DRAG', 'MOMENT_Z');
xlabel('$x/c$');
ylabel('Surface sensitivities');
ylim([-3 2]);
title('Surface sensitivities ITERATION #119')

exportStandardizedFigure(gcf, 'RAE2822_SURFACE_SENSITIVITIES_ITER119', 1, 'legendLocation', 'SouthWest','legendOrientation', 'vertical','addMarkers', false, 'changeColors', false, 'overwriteFigure', true);