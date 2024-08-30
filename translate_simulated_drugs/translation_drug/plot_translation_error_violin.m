clear; clc; close all;

load('all_error_1x_lasso.mat')
load('all_error_2x_lasso.mat')
load('all_error_3x_lasso.mat')
load('all_error_4x_lasso.mat')

%Remove first row since it is the control condition with no drug block
all_error_1x_lasso(1, :) =[];
all_error_2x_lasso(1, :) =[];
all_error_3x_lasso(1, :) =[];
all_error_4x_lasso(1, :) =[];

error_1x = mean(abs(all_error_1x_lasso), 2)*100;
error_2x = mean(abs(all_error_2x_lasso), 2)*100;
error_3x = mean(abs(all_error_3x_lasso), 2)*100;
error_4x = mean(abs(all_error_4x_lasso), 2)*100;
all_error_con = [error_1x; error_2x; error_3x; error_4x];

n_1x = length(error_1x); n_2x = length(error_2x); 
n_3x = length(error_3x); n_4x = length(error_4x);

label_violin = cell((n_1x + n_2x +  n_3x + n_4x), 1); 
label_violin(1:n_1x) = {'1x'}; label_violin(n_1x+1: n_1x +n_2x) = {'2x'};
label_violin(n_1x+n_2x+1: n_1x +n_2x +n_3x) = {'3x'};
label_violin(n_1x+n_2x+n_3x+1: n_1x +n_2x +n_3x +n_4x) = {'4x'};
%color_cell = {[1 0.698 0.4; 1 0.6 0.2; 1 0.5020 0; 0.8 0.4 0]};
color_cell = {[0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5]};

figure(8); hold on; set(gcf, 'color', 'w'); 
vs1 = violinplot(all_error_con, label_violin, 'ViolinColor', color_cell, 'Width', 0.3); ylabel('Translation error %')
%set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
% ylim([0 25]); yticks([ 0 5 10 15 20 25])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(8); set(gcf, 'Units', 'Inches', 'Position', [0 0 6 4], 'PaperUnits', 'Inches', 'PaperSize', [6, 4])
% saveas(gcf,'tranlation_error_avg4.svg')

%% Plot only QT error

all_error_QT = [abs(all_error_1x_lasso(:,2))*100; abs(all_error_2x_lasso(:,2))*100; abs(all_error_3x_lasso(:,2))*100; abs(all_error_4x_lasso(:,2))*100];
color_cell2 = {[0 0 0; 0 0 0; 0 0 0; 0 0 0]};

figure(9); hold on; set(gcf, 'color', 'w'); 
vs2 = violinplot(all_error_QT, label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
ylabel('QTint');
%ylabel('QRSdur'); 
%ylabel('Tpeakend'); 
% ylabel('Twaveamp');
xlabel('ETPC')
ylim([0 25]); yticks([ 0 5 10 15 20 25])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(9); set(gcf, 'Units', 'Inches', 'Position', [0 0 6 4], 'PaperUnits', 'Inches', 'PaperSize', [6, 4])
%f = gcf; exportgraphics(f, 'translation_error.png', 'Resolution', 300)
% saveas(gcf,'QTint_tranlation_error3.svg')

%%