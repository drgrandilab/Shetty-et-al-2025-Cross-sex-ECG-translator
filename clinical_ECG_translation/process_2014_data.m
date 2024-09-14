
clear; close all; clc;
color_m_shade = [1, 0, 1];%[0.968627451 0.305882353 0.839215686]; % Pink
color_f_shade = [0  0.4470  0.7410];%[0.333333333 0.62745098 0.984313725]; % Blue

color_m_solid = [0.8, 0, 0.8];
color_f_solid = [0  0  0.8];

addpath('../build_regression_model/')
load('regression_Blasso4.mat')

% Blasso = diag([1 1 1 1]);

time_point = [-0.5, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 12, 14, 24];

T = readtable('Johannesen_Strauss2014.csv', 'Delimiter',',');
S =  table2struct(T);
Tfilter = table;
variables = T.Properties.VariableNames;
filter_variables = variables(1, [2, 3, 11, 13:24, 26, 30, 33]);

for n = 1: length(filter_variables)
    xvar = T.(string(filter_variables(n)));
    xvar_reshaped = reshape(xvar,[3,length(xvar)/3]);
    if isa(xvar_reshaped,'double')
        xvar_new = mean(xvar_reshaped)';
    else
        xvar_new = xvar_reshaped(1,:);
        xvar_new = xvar_new';
    end
    Tfilter.(string(filter_variables(n))) = xvar_new;
end

%%

ind_dofetilide = find(strcmp(Tfilter.EXTRT, 'Dofetilide'));
ind_quinidine = find(strcmp(Tfilter.EXTRT, 'Quinidine Sulph'));
ind_ranolazine = find(strcmp(Tfilter.EXTRT, 'Ranolazine'));
ind_verapamil = find(strcmp(Tfilter.EXTRT, 'Verapamil HCL'));
ind_placebo = find(strcmp(Tfilter.EXTRT, 'Placebo'));

% Uncomment the desired drug to process
drug_data = Tfilter(ind_dofetilide,:);
%drug_data = Tfilter(ind_quinidine,:);
%drug_data = Tfilter(ind_ranolazine,:);
%drug_data = Tfilter(ind_verapamil,:);

placebo_data = Tfilter(ind_placebo,:);

female_id = 1001:1011;
male_id = 1012:1022;

rel_RR_male = nan(11, 16); rel_qtc_male = nan(11, 16); rel_qrs_male = nan(11, 16); rel_tpeakend_male = nan(11, 16); rel_twaveamp_male = nan(11, 16);
male_master_RR = nan(11,16); male_master_qtc = nan(11, 16); male_master_qrs = nan(11, 16); male_master_tpeakend = nan(11, 16); male_master_twaveamp = nan(11, 16); male_master_drug_conc = nan(11,16);
female_master_RR = nan(11,16); female_master_qtc = nan(11, 16); female_master_qrs = nan(11, 16); female_master_tpeakend = nan(11, 16); female_master_twaveamp = nan(11, 16); female_master_drug_conc = nan(11,16);
translated_master_qtc = nan(11, 16); translated_master_qrs = nan(11, 16); translated_master_tpeakend = nan(11, 16); translated_master_twaveamp = nan(11, 16);


for i = 1:11
    % Male data
    male = male_id(i);
    m_ind = find(drug_data.RANDID== male);
    male_id_RR = table2array(drug_data(m_ind, 12));
    male_id_qtc = table2array(drug_data(m_ind, 18));
    male_id_qrs = table2array(drug_data(m_ind, 15));
    male_id_tpeakend = table2array(drug_data(m_ind, 16));
    male_id_twaveamp = table2array(drug_data(m_ind, 17));
    male_id_drug_conc = table2array(drug_data(m_ind, 10));

    m_ind_placebo = find(placebo_data.RANDID== male);
    male_id_placebo_RR = table2array(placebo_data(m_ind_placebo, 12));
    male_id_placebo_qtc = table2array(placebo_data(m_ind_placebo, 18));
    male_id_placebo_qrs = table2array(placebo_data(m_ind_placebo, 15));
    male_id_placebo_tpeakend = table2array(placebo_data(m_ind_placebo, 16));
    male_id_placebo_twaveamp = table2array(placebo_data(m_ind_placebo, 17));

    %Placebo correction
    rel_RR_male(i,:) = (( (male_id_RR - male_id_RR(1)) - (male_id_placebo_RR - male_id_placebo_RR(1)) )/male_id_RR(1) )*100  ;
    rel_qtc_male(i,:) = (( (male_id_qtc - male_id_qtc(1)) - (male_id_placebo_qtc - male_id_placebo_qtc(1)) )/male_id_qtc(1) )*100  ;
    rel_qrs_male(i,:) = (( (male_id_qrs - male_id_qrs(1)) - (male_id_placebo_qrs - male_id_placebo_qrs(1)) )/male_id_qrs(1) )*100  ;
    rel_tpeakend_male(i,:) = (( (male_id_tpeakend - male_id_tpeakend(1)) - (male_id_placebo_tpeakend - male_id_placebo_tpeakend(1)) )/male_id_tpeakend(1) )*100 ;
    rel_twaveamp_male(i,:) = (( (male_id_twaveamp - male_id_twaveamp(1)) - (male_id_placebo_twaveamp - male_id_placebo_twaveamp(1)) )/male_id_twaveamp(1) )*100  ;

    RR_male_abs = (male_id_RR(1).*(1 + rel_RR_male(i,:)./100));
    qtc_male_abs = (male_id_qtc(1).*(1 + rel_qtc_male(i,:)./100));
    qrs_male_abs = (male_id_qrs(1).*(1 + rel_qrs_male(i,:)./100));
    tpeakend_male_abs = (male_id_tpeakend(1).*(1 + rel_tpeakend_male(i,:)./100));
    twaveamp_male_abs = (male_id_twaveamp(1).*(1 + rel_twaveamp_male(i,:)./100));

    male_master_RR(i, :) = RR_male_abs;
    male_master_qtc(i, :) = qtc_male_abs;
    male_master_qrs(i, :) = qrs_male_abs;
    male_master_tpeakend(i, :) = tpeakend_male_abs;
    male_master_twaveamp(i, :) = twaveamp_male_abs;
    male_master_drug_conc(i,:) = male_id_drug_conc;

end

%%
rel_qrs_male_mean = mean(rel_qrs_male, 'omitnan');
rel_qtc_male_mean = mean(rel_qtc_male, 'omitnan');
rel_tpeakend_male_mean = mean(rel_tpeakend_male, 'omitnan');
rel_twaveamp_male_mean = mean(rel_twaveamp_male, 'omitnan');

% Translated_data
Y_simulated_output = zeros(16, 4);
predicted_percentage_change = zeros(16, 4);
for t = 1:16
    percent_change = [rel_qrs_male_mean(t)  rel_qtc_male_mean(t)   rel_tpeakend_male_mean(t)  rel_twaveamp_male_mean(t) ];
    X_simulated_feature = mean(X_train).*(1+percent_change/100);
    
    Y_baseline_mean = mean(Y_train);
    x = log(X_simulated_feature);
    xz = (x-mean(X_log))./std(X_log);
    yz = xz*Blasso;
    y = yz.*std(Y_log)+mean(Y_log);
    Y_simulated_output(t,:) = exp(y);
    predicted_percentage_change(t,:) =( (Y_simulated_output(t,:) - Y_baseline_mean)./Y_baseline_mean ).*100 ;

end

rel_qtc_translated = predicted_percentage_change(:,2);
rel_qrs_translated = predicted_percentage_change(:,1);
rel_tpeakend_translated = predicted_percentage_change(:,3);
rel_twaveamp_translated = predicted_percentage_change(:,4);

%%  % Female data
    
for j = 1:11
    
    female = female_id(j);
    f_ind = find(drug_data.RANDID== female);
    female_id_RR = table2array(drug_data(f_ind, 12));
    female_id_qtc = table2array(drug_data(f_ind, 18));
    female_id_qrs = table2array(drug_data(f_ind, 15));
    female_id_tpeakend = table2array(drug_data(f_ind, 16));
    female_id_twaveamp = table2array(drug_data(f_ind, 17));
    female_id_drug_conc = table2array(drug_data(f_ind, 10));

    
    f_ind_placebo = find(placebo_data.RANDID== female);
    female_id_placebo_RR = table2array(placebo_data(f_ind_placebo, 12));
    female_id_placebo_qtc = table2array(placebo_data(f_ind_placebo, 18));
    female_id_placebo_qrs = table2array(placebo_data(f_ind_placebo, 15));
    female_id_placebo_tpeakend = table2array(placebo_data(f_ind_placebo, 16));
    female_id_placebo_twaveamp = table2array(placebo_data(f_ind_placebo, 17));
    
    if ~isempty(female_id_qtc)

    %Placebo correction
    rel_RR_female = (( (female_id_RR - female_id_RR(1)) - (female_id_placebo_RR - female_id_placebo_RR(1)) )/female_id_RR(1) )*100  ;    
    rel_qtc_female = (( (female_id_qtc - female_id_qtc(1)) - (female_id_placebo_qtc - female_id_placebo_qtc(1)) )/female_id_qtc(1) )*100  ;
    rel_qrs_female = (( (female_id_qrs - female_id_qrs(1)) - (female_id_placebo_qrs - female_id_placebo_qrs(1)) )/female_id_qrs(1) )*100  ;
    rel_tpeakend_female = (( (female_id_tpeakend - female_id_tpeakend(1)) - (female_id_placebo_tpeakend - female_id_placebo_tpeakend(1)) )/female_id_tpeakend(1) )*100 ;
    rel_twaveamp_female = (( (female_id_twaveamp - female_id_twaveamp(1)) - (female_id_placebo_twaveamp - female_id_placebo_twaveamp(1)) )/female_id_twaveamp(1) )*100  ;

    RR_female_abs = (female_id_RR(1).*(1 + rel_RR_female./100));
    qtc_female_abs = (female_id_qtc(1).*(1 + rel_qtc_female./100));
    qrs_female_abs = (female_id_qrs(1).*(1 + rel_qrs_female./100));
    tpeakend_female_abs = (female_id_tpeakend(1).*(1 + rel_tpeakend_female./100));
    twaveamp_female_abs = (female_id_twaveamp(1).*(1 + rel_twaveamp_female./100));

    female_master_RR(j, :) = RR_female_abs;
    female_master_qtc(j, :) = qtc_female_abs;
    female_master_qrs(j, :) = qrs_female_abs;
    female_master_tpeakend(j, :) = tpeakend_female_abs;
    female_master_twaveamp(j, :) = twaveamp_female_abs;
    female_master_drug_conc(j,:) = female_id_drug_conc;

    %%Translated
    qtc_translated_abs = (female_id_qtc(1).*(1 + rel_qtc_translated./100));
    qrs_translated_abs = (female_id_qrs(1).*(1 + rel_qrs_translated./100));
    tpeakend_translated_abs = (female_id_tpeakend(1).*(1 + rel_tpeakend_translated./100));
    twaveamp_translated_abs = (female_id_twaveamp(1).*(1 + rel_twaveamp_translated./100));

    translated_master_qtc(j, :) = qtc_translated_abs;
    translated_master_qrs(j, :) = qrs_translated_abs;
    translated_master_tpeakend(j, :) = tpeakend_translated_abs;
    translated_master_twaveamp(j, :) = twaveamp_translated_abs;

    end

end


%%

figure(3); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(female_master_qrs, 0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(male_master_qrs, 0.3,color_m_shade,time_point);
[lineOut3, fillOut3] = stdshade(translated_master_qrs, 0.3,'k',time_point);
plot(time_point, mean(female_master_qrs, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(male_master_qrs, 'omitnan'), '.-','Markersize', 12,'linewidth', 1.6, 'Color', color_m_solid)
plot(time_point, mean(translated_master_qrs, 'omitnan'), '.-', 'Markersize', 12, 'linewidth', 1.6, 'Color', 'k')
diff_qrs_translated = abs((mean(female_master_qrs, 'omitnan') - mean(translated_master_qrs, 'omitnan'))./mean(female_master_qrs, 'omitnan'))';
xlim([-1 25]); ylim([75 110]);
xticks([0,6,12,18,24 ]);  yticks([80, 90, 100, 110 ]);
set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
%xlabel('Time (h)'); ylabel('QRS dur (ms)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(3); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2014QRS.svg')
% saveas(gcf,'fig_d1/quin_2014QRS.svg')
% saveas(gcf,'fig_d1/rano_2014QRS.svg')
% saveas(gcf,'fig_d1/vera_2014QRS.svg')

figure(4); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(female_master_qtc,0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(male_master_qtc,0.3,color_m_shade,time_point);
[lineOut3, fillOut3] = stdshade(translated_master_qtc,0.3,'k',time_point);
hold on 
plot(time_point, mean(female_master_qtc, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(male_master_qtc, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_m_solid)
plot(time_point, mean(translated_master_qtc, 'omitnan'), '.-', 'Markersize', 12, 'linewidth', 1.6, 'Color', 'k')
diff_qtc_translated = abs((mean(female_master_qtc, 'omitnan') - mean(translated_master_qtc, 'omitnan'))./mean(female_master_qtc, 'omitnan'))';
xlim([-1 25]); ylim([370 500]); 
xticks([0,6,12,18,24 ]);  yticks([380,420,460,500 ]);
set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
%xlabel('Time (h)'); ylabel('QTc (ms)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(4); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2014QTc.svg')
% saveas(gcf,'fig_d1/quin_2014QTc.svg')
% saveas(gcf,'fig_d1/rano_2014QTc.svg')
% saveas(gcf,'fig_d1/vera_2014QTc.svg')

figure(5); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(female_master_tpeakend,0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(male_master_tpeakend,0.3,color_m_shade,time_point);
[lineOut3, fillOut3] = stdshade(translated_master_tpeakend,0.3,'k',time_point);
plot(time_point, mean(female_master_tpeakend, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(male_master_tpeakend, 'omitnan'), '.-','Markersize', 12,'linewidth', 1.6, 'Color', color_m_solid)
plot(time_point, mean(translated_master_tpeakend, 'omitnan'), '.-', 'Markersize', 12, 'linewidth', 1.6, 'Color', 'k')
diff_tpeakend_translated = abs((mean(female_master_tpeakend, 'omitnan') - mean(translated_master_tpeakend, 'omitnan'))./mean(female_master_tpeakend, 'omitnan'))';
xlim([-1 25]); ylim([60 150]); 
xticks([0,6,12,18,24 ]);  yticks([60,90,120,150]);
set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
%xlabel('Time (h)'); ylabel('T-peak-end dur. (ms)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02],'box', 'off', 'tickdir', 'out');
figure(5); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2014Tpeakend.svg')
% saveas(gcf,'fig_d1/quin_2014Tpeakend.svg')
% saveas(gcf,'fig_d1/rano_2014Tpeakend.svg')
% saveas(gcf,'fig_d1/vera_2014Tpeakend.svg')

figure(6); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(female_master_twaveamp,0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(male_master_twaveamp,0.3,color_m_shade,time_point);
[lineOut3, fillOut3] = stdshade(translated_master_twaveamp,0.3,'k',time_point);
hold on 
plot(time_point, mean(female_master_twaveamp, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(male_master_twaveamp, 'omitnan'), '.-','Markersize', 12,'linewidth', 1.6, 'Color', color_m_solid)
plot(time_point, mean(translated_master_twaveamp, 'omitnan'), '.-', 'Markersize', 12, 'linewidth', 1.6, 'Color', 'k')
diff_twaveamp_translated = abs((mean(female_master_twaveamp, 'omitnan') - mean(translated_master_twaveamp, 'omitnan'))./mean(female_master_twaveamp, 'omitnan'))';
xlim([-1 25]); ylim([250 800]); 
set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
xticks([0,6,12,18,24 ]);  yticks([300,450,600,750]);
%xlabel('Time (h)'); ylabel('T-wave amp. (\muV)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(6); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2014Twaveamp.svg')
% saveas(gcf,'fig_d1/quin_2014Twaveamp.svg')
% saveas(gcf,'fig_d1/rano_2014Twaveamp.svg')
% saveas(gcf,'fig_d1/vera_2014Twaveamp.svg')

female_master_drug_conc(:,1) = 0;   male_master_drug_conc(:,1) = 0;  %Control time point
figure(7); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(female_master_drug_conc,0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(male_master_drug_conc,0.3,color_m_shade,time_point);
hold on 
plot(time_point, mean(female_master_drug_conc, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(male_master_drug_conc, 'omitnan'), '.-','Markersize', 12,'linewidth', 1.6, 'Color', color_m_solid)
xlim([-1 25]); 
set(gca, 'xticklabel', {[]});
xticks([0,6,12,18,24 ]);  
%xlabel('Time (h)'); ylabel('T-wave amp. (\muV)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(7); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2014drug.svg')
% saveas(gcf,'fig_d1/quin_2014drug.svg')
% saveas(gcf,'fig_d1/rano_2014drug.svg')
% saveas(gcf,'fig_d1/vera_2014drug.svg')

figure(8); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(1000./female_master_RR,0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(1000./male_master_RR,0.3,color_m_shade,time_point);
hold on 
plot(time_point, mean(1000./female_master_RR, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(1000./male_master_RR, 'omitnan'), '.-','Markersize', 12,'linewidth', 1.6, 'Color', color_m_solid)
xlim([-1 25]); 
%set(gca, 'xticklabel', {[]});
xticks([0,6,12,18,24 ]);  
%xlabel('Time (h)'); ylabel('T-wave amp. (\muV)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(8); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2014HR.svg')
% saveas(gcf,'fig_d1/quin_2014HR.svg')
% saveas(gcf,'fig_d1/rano_2014HR.svg')
% saveas(gcf,'fig_d1/vera_2014HR.svg')

%%
errors_matrix = [diff_qrs_translated diff_qtc_translated diff_tpeakend_translated diff_twaveamp_translated];
errors_matrix(1, :) = [];  %Ignore first point baseline


