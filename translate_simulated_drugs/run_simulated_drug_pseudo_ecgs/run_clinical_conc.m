% run clinical conc at max plasma difference between male and female
%Quinidine k (male:1820 nM)
% k = [0.3810, 0.8944, 0.8815, 0.9463, 0.6996, 1, 0.80];  %Drug 1

%Quinidine k (female: 2638 nM)
% k = [0.3138, 0.8395, 0.8562, 0.9100, 0.5898, 0.9755, 0.7041]; %Drug 2

%Verapamil k (male: 171.1880 nM)
% k = [0.6272, 0.9762, 0.5451, 1, 0.9705, 1, 1];  %Drug 3

%Verapamil k (female: 334.7304 nM)
k = [0.4625, 0.9546, 0.3643, 1, 0.9505, 1, 1];  %Drug 4

gender = 2;
drug_no = 4;
effect = 'max_conc';

IKr_drug_scale = k(1);
INaL_drug_scale = k(2);
ICaL_drug_scale = k(3);
INa_drug_scale = k(4);
Ito_drug_scale = k(5);
IK1_drug_scale = k(6);
IKs_drug_scale = k(7);

runORd_endo_epi_drug(gender, drug_no, effect, IKr_drug_scale, INaL_drug_scale, ICaL_drug_scale, INa_drug_scale, Ito_drug_scale, IK1_drug_scale, IKs_drug_scale);

