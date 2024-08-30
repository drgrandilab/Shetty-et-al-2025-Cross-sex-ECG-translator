addpath('../')

% Uncomment to plot male baseline cable spacetime
gendertype = 1; num =1;
load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')
load('Vm_gen1_1_CL1000.mat')

% Uncomment to plot female baseline cable spacetime
% gendertype = 2; num =1;
% load('Vm_gen2_1_CL1000.mat')
% load('ts_gen2_Ncell190_endo95_epi95_CL1000.mat')

color_m_solid = [0.8 0 0.8];
color_f_solid = [0  0  0.8];

if gendertype ==1
    idx_cells = 21: 4: 185;
    plot_color = color_m_solid;
elseif gendertype ==2
    idx_cells = 21: 4: 170;
    plot_color = color_f_solid;
end

figure(1); set(gcf, 'color', 'w'); hold on;
for i =  idx_cells 
    plot3(ts_sample - 49000, i.*ones(length(ts_sample), 1), Vm_sample(i, :), 'Color', plot_color, 'Linewidth', 0.5 )
end
% xlabel('Time (ms)');
% zlabel('{\it V_m} (mV)',  "Rotation",0); 
% ylabel('Cell #');
%yticks([1 50 100])
%axis tight
view(3);
%view(0, -45)
xlim([-100 500]);
view(15, -25)
%xlim([-100 700]); ylim([-2 52]); zlim([-90 50])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
%figure(6); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'baseline_space_time_male.svg')
% saveas(gcf,'baseline_space_time_female.svg')
