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
% convert to log 
% all_error_con = log10(all_error_con);

%%

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
%color_cell2 = {[0 0 0; 0 0 0; 0 0 0; 0 0 0]};
color_cell2 = {[0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5]};

figure(9); hold on; set(gcf, 'color', 'w'); 
%vs2 = violinplot(all_error_QT, label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
%vs2 = violinplot(log10(all_error_QT), label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
ylabel('QTint'); yscale log
yticks([0.0001, 0.001, 0.01, 0.1, 1, 10 100]')
%ylabel('QRSdur'); 
%ylabel('Tpeakend'); 
% ylabel('Twaveamp');
xlabel('ETPC')
%ylim([-4 2]); 
ylim([0.0001 100]);
vs2 = violinplot(all_error_QT, label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
% yticks([ 0 5 10 15 20 25])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(9); set(gcf, 'Units', 'Inches', 'Position', [0 0 6 4], 'PaperUnits', 'Inches', 'PaperSize', [6, 4])
%f = gcf; exportgraphics(f, 'translation_error.png', 'Resolution', 300)
% saveas(gcf,'QTint_tranlation_error.svg')
%saveas(gcf,'log_QTint_tranlation_error.svg')
saveas(gcf,'logy_QTint_tranlation_error.svg')
%% Plot only QRS error
all_error_QRS = [abs(all_error_1x_lasso(:,1))*100; abs(all_error_2x_lasso(:,1))*100; abs(all_error_3x_lasso(:,1))*100; abs(all_error_4x_lasso(:,1))*100];
%color_cell2 = {[0 0 0; 0 0 0; 0 0 0; 0 0 0]};
color_cell2 = {[0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5]};

figure(10); hold on; set(gcf, 'color', 'w'); 
%vs3 = violinplot(all_error_QRS, label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
% vs3 = violinplot(log10(all_error_QRS), label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
%ylabel('QTint');
ylabel('QRSdur'); yscale log
yticks([0.0001, 0.001, 0.01, 0.1, 1, 10 100]')
%ylim([-4 2]); 
ylim([0.0001 100]);
vs3 = violinplot(all_error_QRS, label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
%ylabel('Tpeakend'); 
% ylabel('Twaveamp');
xlabel('ETPC')
%ylim([0 25]); yticks([ 0 5 10 15 20 25])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(10); set(gcf, 'Units', 'Inches', 'Position', [0 0 6 4], 'PaperUnits', 'Inches', 'PaperSize', [6, 4])
%f = gcf; exportgraphics(f, 'translation_error.png', 'Resolution', 300)
% saveas(gcf,'QRS_tranlation_error.svg')
%saveas(gcf,'log_QRS_tranlation_error.svg')
saveas(gcf,'logy_QRS_tranlation_error.svg')
%% Plot only T-peak-end error
all_error_tpeakend = [abs(all_error_1x_lasso(:,3))*100; abs(all_error_2x_lasso(:,3))*100; abs(all_error_3x_lasso(:,3))*100; abs(all_error_4x_lasso(:,3))*100];
%color_cell2 = {[0 0 0; 0 0 0; 0 0 0; 0 0 0]};
color_cell2 = {[0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5]};

figure(11); hold on; set(gcf, 'color', 'w'); 
%vs4 = violinplot(all_error_tpeakend, label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
%vs4 = violinplot(log10(all_error_tpeakend), label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
%ylabel('QTint');
%ylabel('QRSdur'); 
ylabel('Tpeakend'); yscale log
yticks([0.0001, 0.001, 0.01, 0.1, 1, 10 100]')
% ylabel('Twaveamp');
xlabel('ETPC');
ylim([0.0001 100]);
vs4 = violinplot(all_error_tpeakend, label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
%ylim([-4 2]); 
%ylim([0 25]); yticks([ 0 5 10 15 20 25])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(11); set(gcf, 'Units', 'Inches', 'Position', [0 0 6 4], 'PaperUnits', 'Inches', 'PaperSize', [6, 4])
%f = gcf; exportgraphics(f, 'translation_error.png', 'Resolution', 300)
% saveas(gcf,'Tpeakend_tranlation_error.svg')
%saveas(gcf,'log_Tpeakend_tranlation_error.svg')
saveas(gcf,'logy_Tpeakend_tranlation_error.svg')
%% Plot only T-wave amp error
all_error_twaveamp = [abs(all_error_1x_lasso(:,4))*100; abs(all_error_2x_lasso(:,4))*100; abs(all_error_3x_lasso(:,4))*100; abs(all_error_4x_lasso(:,4))*100];
%color_cell2 = {[0 0 0; 0 0 0; 0 0 0; 0 0 0]};
color_cell2 = {[0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5]};

figure(12); hold on; set(gcf, 'color', 'w'); 
%vs5 = violinplot(all_error_twaveamp, label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
%vs5 = violinplot(log10(all_error_twaveamp), label_violin, 'ViolinColor', color_cell2, 'Width', 0.3); %ylim([0 4]); 
%ylabel('QTint');
%ylabel('QRSdur'); 
%ylabel('Tpeakend'); 
ylabel('Twaveamp'); yscale log
yticks([0.0001, 0.001, 0.01, 0.1, 1, 10 100]')
xlabel('ETPC');
% ylim([-4 2]); 
ylim([0.0001 100]);
vs5 = violinplot(all_error_twaveamp, label_violin, 'ViolinColor', [0.5, 0.5, 0.5], 'Width', 0.3); %ylim([0 4]); 
%ylim([0 25]); yticks([ 0 5 10 15 20 25])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(12); set(gcf, 'Units', 'Inches', 'Position', [0 0 6 4], 'PaperUnits', 'Inches', 'PaperSize', [6, 4])
%f = gcf; exportgraphics(f, 'translation_error.png', 'Resolution', 300)
% saveas(gcf,'Twaveamp_tranlation_error.svg')
%saveas(gcf,'log_Twaveamp_tranlation_error.svg')
saveas(gcf,'logy_Twaveamp_tranlation_error.svg')

%%