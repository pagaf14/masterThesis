clear;
clc;
close all;

%% Pressure coefficient

surface_flow=readtable('../DIRECT/CFD/surface_flow.csv');

points_ID=surface_flow.PointID;

x_points=surface_flow.x;
x_points_lowerside=x_points(1:100);
x_points_upperside=x_points(101:200);

y_points=surface_flow.y;
y_points_lowerside=y_points(1:100)*3;
y_points_upperside=y_points(101:200)*3;

cp=surface_flow.Pressure_Coefficient;
cp_lowerside=cp(1:100);
cp_upperside=cp(101:200);


figure;
hold on;
plot1=plot(x_points_upperside, y_points_upperside, '-k');
plot2=plot(x_points_lowerside, y_points_lowerside, '-k');

plot3=plot(x_points_upperside, cp_upperside, 'Color', '#0047AB', 'Marker', '*');
plot4=plot(x_points_lowerside, cp_lowerside, 'Color', '#199FD6', 'Marker', '*');

legend([plot3, plot4],'$C_P$ suction side', '$C_P$ pressure side');
xlabel('$x/c$');
ylabel('$C_P$');
set(gca, 'YDir', 'reverse');

exportStandardizedFigure(gcf, 'QuickStart_CP', 1, 'legendLocation', 'NorthEast','legendOrientation', 'vertical','addMarkers', false, 'changeColors', false);

%% Surface sensitivities

surface_adjoint_CONT= readtable('../CONTINUOUS_ADJOINT/CFD/surface_adjoint.csv');
surface_sensitivity_CONT=surface_adjoint_CONT.Surface_Sensitivity;

surface_adjoint_DISC= readtable('../DISCRETE_ADJOINT/CFD/surface_adjoint.csv');
surface_sensitivity_DISC=surface_adjoint_DISC.Surface_Sensitivity;

figure;
hold on;
plot(x_points, surface_sensitivity_CONT, 'Color', '#B100FF', 'Marker', '*');
plot(x_points, surface_sensitivity_DISC, 'Color', '#88ACE0', 'Marker', '*');
legend('SU2 continuous adjoint', 'SU2 discrete adjoint');
xlabel('$x/c$');
ylabel('Surface sensitivities');
ylim([-2 2]);

exportStandardizedFigure(gcf, 'QuickStart_SurfaceSensitivities', 1, 'legendLocation', 'SouthWest','legendOrientation', 'vertical','addMarkers', false, 'changeColors', false);