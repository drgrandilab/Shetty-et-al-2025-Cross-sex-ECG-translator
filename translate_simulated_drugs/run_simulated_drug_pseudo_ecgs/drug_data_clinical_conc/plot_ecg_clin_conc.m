%clear; clc; close all

gendertype = 2; effective = 'max_conc'; num = 4; 

addpath('../../../male_female_cables_pseudo_ecgs/')
load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')
ts= ts_sample;

color_m = [1, 0.4, 1, 0.35];%[0.968627451 0.305882353 0.839215686]; % Pink
color_f = [0.4  0.7  1 0.35] ;%[0.333333333 0.62745098 0.984313725]; % Blue

if gendertype ==1
    Ncell = 165 + 40;
    cell_end = Ncell - 20;
    plot_color = color_m;
elseif gendertype ==2
    Ncell = 150 + 40;
    cell_end = Ncell - 20;
    plot_color = color_f;
end
k = 1; 

ecgname = sprintf('ecg_gen%d_%s_%d.mat',gendertype, effective, num); 
load(ecgname); 

figure(6); hold on; set(gcf, 'color', 'w'); 
plot(ts - 49000, normalized_phi, 'linewidth', 2, 'Color', plot_color); xlim([-10 750]); ylim([-0.3, 1.1])
xlabel('Time (ms)'); ylabel('\phi (a.u)')
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
        
% set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
% figure(6); set(gcf, 'Units', 'Inches', 'Position', [0 0 8 6], 'PaperUnits', 'Inches', 'PaperSize', [8, 6])
