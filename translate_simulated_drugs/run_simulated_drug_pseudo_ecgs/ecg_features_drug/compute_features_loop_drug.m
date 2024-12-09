 clear; clc; close all
addpath('../ecg_data_drug/')

addpath('../../../male_female_cables_pseudo_ecgs/')
load('ts_gen1_Ncell205_endo102_epi103_CL1000.mat')
ts= ts_sample;
effective = '4x';
%% Calculate male features
gendertype = 1;

QRS_dur_4x_male = nan(99,1);
QRS_amp_4x_male = nan(99,1);
QT_int_4x_male = nan(99,1);
ST_avg_4x_male = nan(99,1);
Twave_dur_4x_male = nan(99,1);
T_peakend_dur_4x_male = nan(99,1);
T_wave_amp_4x_male = nan(99,1);
theta_T_4x_male = nan(99,1);

for num = 1:99
    try
        phiname = sprintf('ecg_gen%d_%s_%d.mat', gendertype, effective, num);
        load(phiname); 
        [QRS_dur, QRS_amp, QT_int, ST_avg, Twave_dur, T_peakend_dur, Twave_amp, theta_T] = get_ECG_features(ts, normalized_phi, gendertype, num);
        QRS_dur_4x_male(num,1) = QRS_dur;
        QRS_amp_4x_male(num,1) = QRS_amp;
        QT_int_4x_male(num,1) = QT_int;
        ST_avg_4x_male(num,1) = ST_avg;
        Twave_dur_4x_male(num,1) = Twave_dur;
        T_peakend_dur_4x_male(num,1) = T_peakend_dur;
        T_wave_amp_4x_male(num,1) = Twave_amp;
        theta_T_4x_male(num,1) = theta_T;
    catch
        disp('Not found'); disp(num)

    end
end

QRS_dur_name = sprintf('QRS_dur_4x_male.mat'); save(QRS_dur_name, 'QRS_dur_4x_male');
QRS_amp_name = sprintf('QRS_amp_4x_male.mat'); save(QRS_amp_name, 'QRS_amp_4x_male');
QT_int_name = sprintf('QT_int_4x_male.mat'); save(QT_int_name, 'QT_int_4x_male');
ST_avg_name = sprintf('ST_avg_4x_male.mat'); save(ST_avg_name, 'ST_avg_4x_male');
Twave_dur_name = sprintf('Twave_dur_4x_male.mat'); save(Twave_dur_name, 'Twave_dur_4x_male');
T_peakend_dur_name = sprintf('T_peakend_dur_4x_male.mat'); save(T_peakend_dur_name, 'T_peakend_dur_4x_male');
T_wave_amp_name = sprintf('T_wave_amp_4x_male.mat'); save(T_wave_amp_name, 'T_wave_amp_4x_male');
theta_T_name = sprintf('theta_T_4x_male.mat'); save(theta_T_name, 'theta_T_4x_male');

%% Calculate female features
gendertype = 2;

QRS_dur_4x_female = nan(99,1);
QRS_amp_4x_female = nan(99,1);
QT_int_4x_female = nan(99,1);
ST_avg_4x_female = nan(99,1);
Twave_dur_4x_female = nan(99,1);
T_peakend_dur_4x_female = nan(99,1);
T_wave_amp_4x_female = nan(99,1);
theta_T_4x_female = nan(99,1);

for num = 1:99
    try
        phiname = sprintf('ecg_gen%d_%s_%d.mat', gendertype, effective, num);
        load(phiname); 
        [QRS_dur, QRS_amp, QT_int, ST_avg, Twave_dur, T_peakend_dur, Twave_amp, theta_T] = get_ECG_features(ts, normalized_phi, gendertype, num);
        QRS_dur_4x_female(num,1) = QRS_dur;
        QRS_amp_4x_female(num,1) = QRS_amp;
        QT_int_4x_female(num,1) = QT_int;
        ST_avg_4x_female(num,1) = ST_avg;
        Twave_dur_4x_female(num,1) = Twave_dur;
        T_peakend_dur_4x_female(num,1) = T_peakend_dur;
        T_wave_amp_4x_female(num,1) = Twave_amp;
        theta_T_4x_female(num,1) = theta_T;
    catch
        disp('Not found'); disp(num)

    end
end

QRS_dur_name = sprintf('QRS_dur_4x_female.mat'); save(QRS_dur_name, 'QRS_dur_4x_female');
QRS_amp_name = sprintf('QRS_amp_4x_female.mat'); save(QRS_amp_name, 'QRS_amp_4x_female');
QT_int_name = sprintf('QT_int_4x_female.mat'); save(QT_int_name, 'QT_int_4x_female');
ST_avg_name = sprintf('ST_avg_4x_female.mat'); save(ST_avg_name, 'ST_avg_4x_female');
Twave_dur_name = sprintf('Twave_dur_4x_female.mat'); save(Twave_dur_name, 'Twave_dur_4x_female');
T_peakend_dur_name = sprintf('T_peakend_dur_4x_female.mat'); save(T_peakend_dur_name, 'T_peakend_dur_4x_female');
T_wave_amp_name = sprintf('T_wave_amp_4x_female.mat'); save(T_wave_amp_name, 'T_wave_amp_4x_female');
theta_T_name = sprintf('theta_T_4x_female.mat'); save(theta_T_name, 'theta_T_4x_female');