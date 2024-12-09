clear; clc; close all
T = readtable('errors_avg_clinical_data.csv', 'Delimiter',',');

addpath('../translate_simulated_drugs/translation_drug/')
diff_qtc = table2array(T(:, 3));
diff_qrs = table2array(T(:, 2));
diff_tpeakend = table2array(T(:, 4));
diff_twaveamp = table2array(T(:, 5));
% idx_drug = 1:120;
% idx_drug = 2:16;  %dof2014
% idx_drug = 18:32; %quin2014
% idx_drug = 34:48; %rano2014
% idx_drug = 50:64; %vera2014
% idx_drug = 66:78;  %dof2016
% idx_drug = 80:92;  %lid_dof2016
% idx_drug = 94:106; %mox_dil2016
% idx_drug = 108:120;  %mex_dof2016

%%
diff_qtc = diff_qtc(idx_drug,:)*100;
diff_qrs = diff_qrs(idx_drug,:)*100;
diff_tpeakend = diff_tpeakend(idx_drug,:)*100;
diff_twaveamp = diff_twaveamp(idx_drug,:)*100;

diff_qtc = diff_qtc(~isnan(diff_qtc));
diff_qrs = diff_qrs(~isnan(diff_qrs));
diff_tpeakend = diff_tpeakend(~isnan(diff_tpeakend));
diff_twaveamp = diff_twaveamp(~isnan(diff_twaveamp));

all_error_con = [diff_qrs; diff_qtc; diff_tpeakend; diff_twaveamp];

n_1x = length(idx_drug); n_2x = length(idx_drug); 
n_3x = length(idx_drug); n_4x = length(idx_drug);

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
yscale log
yticks([0.0001, 0.001, 0.01, 0.1, 1, 10 100]')
ylim([0.0001 100]);
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(8); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'dof2014.svg')
% saveas(gcf,'quin2014.svg')
% saveas(gcf,'rano2014.svg')
% saveas(gcf,'vera2014.svg')
% saveas(gcf,'dof2016.svg')
% saveas(gcf,'lid_dof2016.svg')
% saveas(gcf,'mox_dil2016.svg')
% saveas(gcf,'mex_dof2016.svg')
