clear;
clc;
close all;

filename = 'ONERAM6_HISTORY.pdf'; % specifica il nome del file che desideri cancellare
fullFilePath = fullfile(pwd, filename); % costruisci il percorso completo del file

if exist(fullFilePath, 'file') == 2 % verifica se il file esiste
    delete(fullFilePath); % cancella il file
    disp(['Il file ', filename, ' è stato cancellato con successo.']);
else
    disp(['Il file ', filename, ' non esiste nella cartella corrente.']);
end

% OSSERVAZIONI
% Non so come plottare gli FFD box;
% Il target di C_L=0.286 è rispettato ma il drag non si è ancora assestato
% su un plateu, aumentare il numero di iterazioni;

%% Drag and lift history

history_project=readtable("../history_project.csv");
history_eval=history_project.EVALUATION;
history_drag=history_project.DRAG;
history_lift=history_project.LIFT;

figure;
hold on;

%yyaxis left; % Utilizza l'asse y sinistro
plot(history_eval, history_drag, 'Color', '#C4151C', 'Marker', '*', 'LineWidth', 1.5);
ylabel('$C_D$');
ylim([min(history_drag), max(history_drag)]);
set(gca, 'YColor', 'k');

% yyaxis right; % Utilizza l'asse y destro
% plot(history_eval, history_lift, 'Color', '#0047AB', 'Marker', '*', 'LineWidth', 1.5);
% ylabel('$C_L$');
% ylim([min(history_lift), max(history_lift)]);
% set(gca, 'YColor', 'k');

xlabel('Major Optimizer Iterations');
%legend('$C_D$');%, '$C_L$')

exportStandardizedFigure(gcf, 'ONERAM6_HISTORY', 1, 'legendLocation', 'NorthEast','legendOrientation', 'vertical', 'addMarkers', false, 'changeColors', false, 'changeLineStyle', false);
% 
% close;
% 
% filename = 'ONERAM6_HISTORY.pdf'; % specifica il nome del file PDF che desideri aprire
% fullFilePath = fullfile(pwd, filename); % costruisci il percorso completo del file
% 
% if exist(fullFilePath, 'file') == 2 % verifica se il file esiste
%     open(fullFilePath); % apri il file PDF
% else
%     disp(['Il file ', filename, ' non esiste nella cartella corrente.']);
% end