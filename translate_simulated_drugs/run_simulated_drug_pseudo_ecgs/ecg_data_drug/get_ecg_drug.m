%%
addpath('../../../male_female_cables_pseudo_ecgs/')


gendertype = 1; effective = '4x';
CV_all_4x_male = nan(99,1);
load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')
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
k = 1;
for num = 1:99
    try
        Vname = sprintf('Vm_gen%d_%s_%d.mat',gendertype, effective, num);
        CVname = sprintf('cv_cable_gen%d_%s_%d.mat',gendertype, effective, num);
        load(Vname); 
        load(CVname);
        % 
        ts= ts_sample; Vm = Vm_sample;
        phi = zeros(length(ts), 1);
     
        for t = 1: length(ts)
            for cellno = 21:cell_end
                gradV = Vm(cellno - 1, t) - Vm(cellno + 1, t);
                r = 2 + (Ncell*0.01 -cellno*0.01);
                phi(t) = phi(t)+ (k*gradV/(r*r))*0.01 ;
            end
        end
        normalized_phi = phi/0.4446;

   
        CV_all_4x_male(num) = cv_est(end);
        figure(6); hold on; set(gcf, 'color', 'w'); 
        plot(ts - 49000, normalized_phi, 'linewidth', 1.5, 'Color', plot_color); xlim([-10 700])
        xlabel('Time (ms)'); ylabel('phi (a.u)')
        set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
        % ecgname = sprintf('ecg_gen%d_%s_%d.mat',gendertype, effective, num); save(ecgname, 'normalized_phi');

     catch
        disp('Not found'); disp(num)
     end
end

%%
% figure(6); set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 6], 'PaperUnits', 'Inches', 'PaperSize', [8, 6])
% f = gcf; exportgraphics(f,'female_population.png','Resolution', 300)

