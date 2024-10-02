clear; clc; close all;

output_names = {'CV', 'QRSamp', 'QRSdur', 'QTint', 'STavg',...
    'Tpeakend dur', 'Tawave amp', 'theta Twave', 'Twave dur'};

load('weird_ecg_Jan.mat');

load('CV_all_male.mat'); [CV_all_male,~] = removerows(CV_all_male,'ind', weird_ecg_Jan);
load('QRS_amp_male.mat'); [QRS_amp_male,~] = removerows(QRS_amp_male,'ind', weird_ecg_Jan);
load('QRS_dur_male.mat'); [QRS_dur_male,~] = removerows(QRS_dur_male,'ind', weird_ecg_Jan);
load('QT_int_male.mat'); [QT_int_male,~] = removerows(QT_int_male,'ind', weird_ecg_Jan);
load('ST_avg_male.mat'); [ST_avg_male,~] = removerows(ST_avg_male,'ind', weird_ecg_Jan);
load('T_peakend_dur_male.mat'); [T_peakend_dur_male,~] = removerows(T_peakend_dur_male,'ind', weird_ecg_Jan);
load('T_wave_amp_male.mat'); [T_wave_amp_male,~] = removerows(T_wave_amp_male,'ind', weird_ecg_Jan);
load('theta_T_male.mat'); [theta_T_male,~] = removerows(theta_T_male,'ind', weird_ecg_Jan);
load('Twave_dur_male.mat'); [Twave_dur_male,~] = removerows(Twave_dur_male,'ind', weird_ecg_Jan);

load('CV_all_female.mat'); [CV_all_female,~] = removerows(CV_all_female,'ind', weird_ecg_Jan);
load('QRS_amp_female.mat'); [QRS_amp_female,~] = removerows(QRS_amp_female,'ind', weird_ecg_Jan);
load('QRS_dur_female.mat'); [QRS_dur_female,~] = removerows(QRS_dur_female,'ind', weird_ecg_Jan);
load('QT_int_female.mat'); [QT_int_female,~] = removerows(QT_int_female,'ind', weird_ecg_Jan);
load('ST_avg_female.mat'); [ST_avg_female,~] = removerows(ST_avg_female,'ind', weird_ecg_Jan);
load('T_peakend_dur_female.mat'); [T_peakend_dur_female,~] = removerows(T_peakend_dur_female,'ind', weird_ecg_Jan);
load('T_wave_amp_female.mat'); [T_wave_amp_female,~] = removerows(T_wave_amp_female,'ind', weird_ecg_Jan);
load('theta_T_female.mat'); [theta_T_female,~] = removerows(theta_T_female,'ind', weird_ecg_Jan);
load('Twave_dur_female.mat'); [Twave_dur_female,~] = removerows(Twave_dur_female,'ind', weird_ecg_Jan);


CV = [CV_all_female; CV_all_male];
QRSamp = [QRS_amp_female; QRS_amp_male];
QRS_dur = [QRS_dur_female; QRS_dur_male];
QT_interval = [QT_int_female; QT_int_male];
ST_avg = [ST_avg_female; ST_avg_male];
Tpeak_end = [T_peakend_dur_female; T_peakend_dur_male];
Twave_amp = [T_wave_amp_female; T_wave_amp_male];
theta_Twave = [theta_T_female; theta_T_male];
Twave_dur = [Twave_dur_female; Twave_dur_male];

n_samples = length(CV_all_male);
label_violin = cell(2*n_samples, 1); label_violin(1:n_samples) = {'Female'}; label_violin(n_samples +1:end) = {'Male'};
color_cell = {[0.333333333 0.62745098 0.984313725; 0.968627451 0.305882353 0.839215686]};

figure(8); hold on; set(gcf, 'color', 'w'); 
subplot(3,3,1); hold on; plot(nan, nan, '.','color', [1 1 1]); vs1 = violinplot(CV , label_violin, 'ViolinColor', color_cell, 'Width', 0.3); legend('CV (cm/s)'); legend boxoff 
subplot(3,3,2); hold on; plot(nan, nan, '.', 'color', [1 1 1]); vs2 = violinplot(QRSamp , label_violin, 'ViolinColor', color_cell, 'Width', 0.3); legend('QRS amp. (a.u)'); legend boxoff 
subplot(3,3,3); hold on; plot(nan, nan, '.','color', [1 1 1]); vs3 = violinplot(QRS_dur , label_violin, 'ViolinColor', color_cell, 'Width', 0.3); legend('QRS dur. (ms)'); legend boxoff 
subplot(3,3,4); hold on; plot(nan, nan,'.', 'color', [1 1 1]); vs4 = violinplot(QT_interval , label_violin, 'ViolinColor', color_cell, 'Width', 0.3); legend('QT interval (ms)'); legend boxoff 
subplot(3,3,5); hold on; plot(nan, nan,'.', 'color', [1 1 1]); vs5 = violinplot(ST_avg , label_violin, 'ViolinColor', color_cell, 'Width', 0.3); legend('|Avg ST_{vol}(a.u)|'); legend boxoff 
subplot(3,3,6); hold on; plot(nan, nan, '.','color', [1 1 1]); vs6 = violinplot(Tpeak_end , label_violin, 'ViolinColor', color_cell, 'Width', 0.3); legend('T-peak-end dur. (ms)'); legend boxoff 
subplot(3,3,7); hold on; plot(nan, nan, '.','color', [1 1 1]); vs7 = violinplot(Twave_amp , label_violin, 'ViolinColor', color_cell, 'Width', 0.3); legend('T-wave amp. (a.u)'); legend boxoff 
subplot(3,3,8); hold on; plot(nan, nan,'.', 'color', [1 1 1]); vs8 = violinplot(theta_Twave , label_violin, 'ViolinColor', color_cell, 'Width', 0.3);legend('T-wave angle (\circ)'); legend boxoff 
subplot(3,3,9); hold on; plot(nan, nan, '.','color', [1 1 1]); vs9 = violinplot(Twave_dur , label_violin, 'ViolinColor', color_cell, 'Width', 0.3); legend('T-wave dur. (ms)'); legend boxoff 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(8); set(gcf, 'Units', 'Inches', 'Position', [0 0 12 10], 'PaperUnits', 'Inches', 'PaperSize', [12, 10])
    %f = gcf; exportgraphics(f, 'feature_distribution.png', 'Resolution', 300)

