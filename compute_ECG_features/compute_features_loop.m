clear; clc; close all
addpath('../')
addpath('/media/roshni/DATADRIVE01/Roshni/Cross-sex-Translator-RS/male_female_cables_pseudo_ecgs/plot_pseudo_ecg_populations/ecg_data/');

load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')
ts = ts_sample;

%% Compute male features
gendertype = 1; 
QRS_dur_male = nan(1000,1);
QRS_amp_male = nan(1000,1);
QT_int_male = nan(1000,1);
ST_avg_male = nan(1000,1);
Twave_dur_male = nan(1000,1);
T_peakend_dur_male = nan(1000,1);
T_wave_amp_male = nan(1000,1);
theta_T_male = nan(1000,1);

for num = 1:1000
    try
        phiname = sprintf('ecg_gen%d_%d.mat', gendertype, num);
        load(phiname); 
        [QRS_dur, QRS_amp, QT_int, ST_avg, Twave_dur, T_peakend_dur, Twave_amp, theta_T] = get_ECG_features(ts, normalized_phi, gendertype, num);
        QRS_dur_male(num,1) = QRS_dur;
        QRS_amp_male(num,1) = QRS_amp;
        QT_int_male(num,1) = QT_int;
        ST_avg_male(num,1) = ST_avg;
        Twave_dur_male(num,1) = Twave_dur;
        T_peakend_dur_male(num,1) = T_peakend_dur;
        T_wave_amp_male(num,1) = Twave_amp;
        theta_T_male(num,1) = theta_T;
    catch
        disp('Not found'); disp(num)
         
    end
end


QRS_dur_name = sprintf('QRS_dur_male.mat'); save(QRS_dur_name, 'QRS_dur_male');
QRS_amp_name = sprintf('QRS_amp_male.mat'); save(QRS_amp_name, 'QRS_amp_male');
QT_int_name = sprintf('QT_int_male.mat'); save(QT_int_name, 'QT_int_male');
ST_avg_name = sprintf('ST_avg_male.mat'); save(ST_avg_name, 'ST_avg_male');
Twave_dur_name = sprintf('Twave_dur_male.mat'); save(Twave_dur_name, 'Twave_dur_male');
T_peakend_dur_name = sprintf('T_peakend_dur_male.mat'); save(T_peakend_dur_name, 'T_peakend_dur_male');
T_wave_amp_name = sprintf('T_wave_amp_male.mat'); save(T_wave_amp_name, 'T_wave_amp_male');
theta_T_name = sprintf('theta_T_male.mat'); save(theta_T_name, 'theta_T_male');

%% Compute female features
gendertype = 2; 
QRS_dur_female = nan(1000,1);
QRS_amp_female = nan(1000,1);
QT_int_female = nan(1000,1);
ST_avg_female = nan(1000,1);
Twave_dur_female = nan(1000,1);
T_peakend_dur_female = nan(1000,1);
T_wave_amp_female = nan(1000,1);
theta_T_female = nan(1000,1);

for num = 1:1000
    try
        phiname = sprintf('ecg_gen%d_%d.mat', gendertype, num);
        load(phiname); 
        [QRS_dur, QRS_amp, QT_int, ST_avg, Twave_dur, T_peakend_dur, Twave_amp, theta_T] = get_ECG_features(ts, normalized_phi, gendertype, num);
        QRS_dur_female(num,1) = QRS_dur;
        QRS_amp_female(num,1) = QRS_amp;
        QT_int_female(num,1) = QT_int;
        ST_avg_female(num,1) = ST_avg;
        Twave_dur_female(num,1) = Twave_dur;
        T_peakend_dur_female(num,1) = T_peakend_dur;
        T_wave_amp_female(num,1) = Twave_amp;
        theta_T_female(num,1) = theta_T;
    catch
        disp('Not found'); disp(num)
         
    end
end


QRS_dur_name = sprintf('QRS_dur_female.mat'); save(QRS_dur_name, 'QRS_dur_female');
QRS_amp_name = sprintf('QRS_amp_female.mat'); save(QRS_amp_name, 'QRS_amp_female');
QT_int_name = sprintf('QT_int_female.mat'); save(QT_int_name, 'QT_int_female');
ST_avg_name = sprintf('ST_avg_female.mat'); save(ST_avg_name, 'ST_avg_female');
Twave_dur_name = sprintf('Twave_dur_female.mat'); save(Twave_dur_name, 'Twave_dur_female');
T_peakend_dur_name = sprintf('T_peakend_dur_female.mat'); save(T_peakend_dur_name, 'T_peakend_dur_female');
T_wave_amp_name = sprintf('T_wave_amp_female.mat'); save(T_wave_amp_name, 'T_wave_amp_female');
theta_T_name = sprintf('theta_T_female.mat'); save(theta_T_name, 'theta_T_female');
