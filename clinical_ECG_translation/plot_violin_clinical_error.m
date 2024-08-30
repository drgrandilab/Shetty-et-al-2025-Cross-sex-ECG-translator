clear; clc; close all
T = readtable('errors_avg_clinical_data.csv', 'Delimiter',',');

addpath('../translate_simulated_drugs/translation_drug/')
diff_qtc = table2array(T(:, 3));
diff_qrs = table2array(T(:, 2));
diff_tpeakend = table2array(T(:, 4));
diff_twaveamp = table2array(T(:, 5));

diff_qtc = diff_qtc(~isnan(diff_qtc))*100;
diff_qrs = diff_qrs(~isnan(diff_qrs))*100;
diff_tpeakend = diff_tpeakend(~isnan(diff_tpeakend))*100;
diff_twaveamp = diff_twaveamp(~isnan(diff_twaveamp))*100;

all_error_con = [diff_qrs; diff_qtc; diff_tpeakend; diff_twaveamp];

n_1x = 112; n_2x = 112; 
n_3x = 112; n_4x = 112;

label_violin = cell((n_1x + n_2x +  n_3x + n_4x), 1); 
label_violin(1:n_1x) = {'QRSdur'}; label_violin(n_1x+1: n_1x +n_2x) = {'QTc'};
label_violin(n_1x+n_2x+1: n_1x +n_2x +n_3x) = {'Tpeakend'};
label_violin(n_1x+n_2x+n_3x+1: n_1x +n_2x +n_3x +n_4x) = {'Twaveamp'};
%color_cell = {[1 0.698 0.4; 1 0.6 0.2; 1 0.5020 0; 0.8 0.4 0]};
color_cell = {[0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5]};

figure(8); hold on; set(gcf, 'color', 'w'); 
vs1 = violinplot(all_error_con, label_violin, 'ViolinColor', color_cell, 'Width', 0.3); ylabel('Translation error %')
%set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
%ylim([0 13]); yticks([ 0 3 6 9 12])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1.2, 'box', 'off', 'tickdir', 'out');
figure(8); set(gcf, 'Units', 'Inches', 'Position', [0 0 6 4], 'PaperUnits', 'Inches', 'PaperSize', [6, 4])
saveas(gcf,'clinical_error_violin.svg')