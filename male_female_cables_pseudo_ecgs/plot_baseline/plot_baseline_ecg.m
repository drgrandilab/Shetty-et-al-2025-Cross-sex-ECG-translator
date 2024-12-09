addpath('../')

% Uncomment to plot baseline male ecg
gendertype = 2; num =1;
load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')
load('Vm_gen1_1_CL1000.mat')

% Uncomment to plot baseline female ecg
% gendertype = 2; num =1;
% load('Vm_gen2_1_CL1000.mat')
% load('ts_gen2_Ncell190_endo95_epi95_CL1000.mat')


ts = ts_sample;
Vm = Vm_sample;

phi = zeros(length(ts), 1);
k = 1;
color_m = [0.968627451 0.305882353 0.839215686]; % Pink
color_f = [0.333333333 0.62745098 0.984313725]; % Blue

if gendertype ==1
    Ncell = 165 + 40;
    cell_end = Ncell - 20;
    plot_color = color_m;
elseif gendertype ==2
    Ncell = 150 + 40;
    cell_end = Ncell - 20;
    plot_color = color_f;
end

for t = 1: length(ts)
    for cellno = 21:cell_end
        gradV = Vm(cellno - 1, t) - Vm(cellno + 1, t);
        r = 2 + (Ncell*0.01 -cellno*0.01);
        phi(t) = phi(t)+ (k*gradV/(r*r))*0.01 ;
    end
end
normalized_phi = phi/0.4446;   %Normailize to 0.4446; peak value for baseline male


figure(6); hold on; set(gcf, 'color', 'w'); 
plot(ts - 49000, normalized_phi, 'linewidth', 1.5, 'Color', plot_color); xlim([-10 500])
xlabel('Time (ms)'); ylabel('phi (a.u)')
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');

% ecgname = sprintf('ecg_gen_baseline_gen%d.mat',gendertype); save(ecgname, 'normalized_phi');
