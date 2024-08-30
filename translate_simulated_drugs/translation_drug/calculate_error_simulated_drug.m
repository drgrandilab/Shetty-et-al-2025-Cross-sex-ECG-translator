clear; close all; clc
color = [0 0 0]; % Black
color_m = [0.968627451 0.305882353 0.839215686]; % Pink
color_f = [0.333333333 0.62745098 0.984313725]; % Blue

addpath('../run_simulated_drug_pseudo_ecgs/ecg_features_drug/')

failed_1x = [39 56 65 71 94 95];
failed_2x = [2 14 25 34  39  45  52  56  65  68 71  72  73  94  95  97 ];
failed_3x = [2 14 25 34  39  42  45  46 52  56 62 65  68 71  72  73  94  95 96 97 98];
failed_4x = [2 8 14 21 24 25 34  39 41  42  45  46 52  56  61  62 65 67 68 71  72  73  94  95 96 97 98];

load('Y_pred_4x_lasso.mat')

load('QRS_dur_4x_female.mat'); [QRS_dur_4x_female,~] = removerows(QRS_dur_4x_female,'ind', failed_4x);
load('QT_int_4x_female.mat'); [QT_int_4x_female,~] = removerows(QT_int_4x_female,'ind', failed_4x);
load('T_peakend_dur_4x_female.mat'); [T_peakend_dur_4x_female,~] = removerows(T_peakend_dur_4x_female,'ind', failed_4x);
load('T_wave_amp_4x_female.mat'); [T_wave_amp_4x_female,~] = removerows(T_wave_amp_4x_female,'ind', failed_4x);

Y_true_drug= [QRS_dur_4x_female, QT_int_4x_female, T_peakend_dur_4x_female, T_wave_amp_4x_female];

output_names = {'QRSdur', 'QTint', 'Tpeakend dur', 'Tawave amp'};

test_count = length(QRS_dur_4x_female);
N_outputs_Y = length(output_names);

all_error_4x_lasso = zeros(test_count, N_outputs_Y);
% (1 − predicted value/actual value)
for feature =1:4
    all_error_4x_lasso(:, feature) = 1 - Y_pred_4x_lasso(:, feature)./Y_true_drug(:, feature);

end