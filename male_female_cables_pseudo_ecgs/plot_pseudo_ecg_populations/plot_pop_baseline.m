clear; close all; clc;
addpath('../')

color_m_shade = [1, 0, 1, 0.03];    % Pink
color_f_shade = [0  0  0.8, 0.03];  % Blue

color_m_solid = [1, 0, 1  1];       % Pink
color_f_solid = [0  0  0.8  1];     % Blue

time =  cell2mat(struct2cell(load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')));
% gendertype = 1;
gendertype = 2;

ts = time - 49000;
if gendertype ==1
    plot_color = color_m_solid;
    plot_color2 = color_m_shade;
    baseline_ECG = cell2mat(struct2cell(load('ecg_gen_baseline_gen1.mat')));

elseif gendertype ==2
    plot_color = color_f_solid;
    plot_color2 = color_f_shade;
    baseline_ECG = cell2mat(struct2cell(load('ecg_gen_baseline_gen2.mat')));
end

figure(1); hold on; set(gcf, 'color', 'w'); 

for k = 1:400   %Plotting first 400
    try
    ecg_name_i = sprintf('ecg_data/ecg_gen%d_%d.mat',  gendertype, k);    %Plot only cables wih propagating APs
    ecg_plot_i = cell2mat(struct2cell(load(ecg_name_i)));
    plot(ts, ecg_plot_i, 'linewidth', 0.5, 'Color', plot_color2); xlim([-10 600]); ylim([-0.1, 1.1])

    catch
        disp(k)
    end

end

figure(1); hold on; 
plot(ts, baseline_ECG, 'linewidth', 3, 'Color', plot_color); xlim([-10 600]); ylim([-0.1, 1.1])
xlabel('Time (ms)'); ylabel('\phi (a.u)')
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
% saveas(gcf,'pop_baseline_pseudoECG_male.svg')
% saveas(gcf,'pop_baseline_pseudoECG_female.svg')
