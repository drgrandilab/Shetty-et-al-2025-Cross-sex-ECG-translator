%Dof CIPA 1x   %ECGNo.30  

addpath('../../male_female_cables_pseudo_ecgs/')
load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')

addpath('../run_simulated_drug_pseudo_ecgs/ecg_data_drug/')
male_control =  cell2mat(struct2cell(load('ecg_gen1_1x_1.mat')));
male_dof =  cell2mat(struct2cell(load('ecg_gen1_1x_30.mat')));
female_dof =  cell2mat(struct2cell(load('ecg_gen2_1x_30.mat')));

figure(1); hold on; set(gcf, 'color', 'w'); 
plot(ts_sample - 49000, male_control , 'linewidth', 2, 'Color', [1, 0, 1]); 
plot(ts_sample - 49000, male_dof , 'linewidth', 2, 'Color', [0.6431, 0.2549, 0.7020]); 
plot(ts_sample - 49000, female_dof , 'linewidth', 2, 'Color', [0, 0, 0.8]); xlim([-10 650])
xlabel('Time (ms)'); ylabel('phi (a.u)')
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');


%54.2773494649994	410.783496524407	41.0059250786227	0.121854579070435
qrs_dur = 54.2773494649994;
qt_int = 410.783496524407;
tpeak_end_dur = 41.0059250786227;
Twave_amp = 0.121854579070435;

plot(0, 0, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
%plot(36.03, qrs_amp, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
plot(qrs_dur, 0.012, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
plot(qt_int, 0, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
plot(qt_int - tpeak_end_dur, Twave_amp, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
%plot(qt_int - Twave_dur, 0, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
%xlabel('Time (ms)'); ylabel('phi (a.u)')
saveas(gcf,'dofetilide_pred.svg')