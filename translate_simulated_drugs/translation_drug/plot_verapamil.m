%Vera CIPA 1x   %ECGNo.98 (after deletion of failed no.92) 

addpath('../../male_female_cables_pseudo_ecgs/')
load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')

addpath('../run_simulated_drug_pseudo_ecgs/ecg_data_drug/')
male_control =  cell2mat(struct2cell(load('ecg_gen1_1x_1.mat')));
male_vera =  cell2mat(struct2cell(load('ecg_gen1_1x_98.mat')));
female_vera =  cell2mat(struct2cell(load('ecg_gen2_1x_98.mat')));

figure(1); hold on; set(gcf, 'color', 'w'); 
plot(ts_sample - 49000, male_control , 'linewidth', 2, 'Color', [1, 0, 1]); 
plot(ts_sample - 49000, male_vera , 'linewidth', 2, 'Color', [0.6431, 0.2549, 0.7020]); 
plot(ts_sample - 49000, female_vera , 'linewidth', 2, 'Color', [0, 0, 0.8]); xlim([-10 650])
xlabel('Time (ms)'); ylabel('phi (a.u)')
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');


%data_pred = [54.2773494649994	410.783496524407	41.0059250786227    0.121854579070435]  %Dof
%data_pred = [56.6646358878144	607.224982706314	26.6981962195465	0.0666654989714234];%Quin
%data_pred = [54.1632668262273	380.226984436851	40.1538111979159	0.121889033856542]; %Rano
data_pred =  [51.0631531947162	377.900494813694	35.9116138349004	0.102944187584645]; %Vera
qrs_dur = data_pred(1);
qt_int = data_pred(2);
tpeak_end_dur = data_pred(3);
Twave_amp = data_pred(4);

plot(0, 0, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
%plot(36.03, qrs_amp, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
plot(qrs_dur, 0.012, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
plot(qt_int, 0, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
plot(qt_int - tpeak_end_dur, Twave_amp, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
%plot(qt_int - Twave_dur, 0, 'o', 'MarkerSize',12, 'MarkerEdgeColor', 'k')
%xlabel('Time (ms)'); ylabel('phi (a.u)')
saveas(gcf,'vera_pred.svg')