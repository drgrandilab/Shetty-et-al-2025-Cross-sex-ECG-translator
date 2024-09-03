clear; close all; clc;

color_m_shade = [1, 0, 1];%[0.968627451 0.305882353 0.839215686]; % Pink
color_f_shade = [0  0.4470  0.7410];%[0.333333333 0.62745098 0.984313725]; % Blue

color_m_solid = [0.8, 0, 0.8];
color_f_solid = [0  0  0.8];


addpath('../build_regression_model/')
load('regression_Blasso4.mat')

time_point = [-0.5, 1.5, 2, 2.5, 3, 6.5, 7, 7.5, 8, 12, 12.5, 13, 13.5, 24];

T = readtable('Johanessen_Strauss2016.csv', 'Delimiter',',');
S =  table2struct(T);
Tfilter = table;
variables = T.Properties.VariableNames;
filter_variables = variables(1, [2, 3, 11, 13:28, 31, 34]);

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
ind_dofetilide = find(strcmp(Tfilter.TRTA, 'Dofetilide'));
ind_lid_dof = find(strcmp(Tfilter.TRTA, 'Lidocaine + Dofetilide'));
ind_mox_dilt = find(strcmp(Tfilter.TRTA, 'Moxifloxacin + Diltiazem'));
ind_mex_dof = find(strcmp(Tfilter.TRTA, 'Mexiletine + Dofetilide'));
ind_placebo = find(strcmp(Tfilter.TRTA, 'Placebo'));

% Uncomment the desired drug to process
drug_data = Tfilter(ind_dofetilide,:);
% drug_data = Tfilter(ind_lid_dof,:);
% drug_data = Tfilter(ind_mox_dilt,:);
% drug_data = Tfilter(ind_mex_dof,:);
placebo_data = Tfilter(ind_placebo,:);

%To plot drug concentrations
%DOF:5; LIDO:6; MEXI: 7; MOXI:8; MOXI.M2: 9; DILT: 10;
drug1_c = 5; drug2_c = 5 ; %Uncomment for ind_dofetilide
% drug1_c = 5; drug2_c = 6 ; %Uncomment for ind_lid_dof
% drug1_c = 8; drug2_c = 10 ; %Uncomment for ind_mox_dilt
% drug1_c = 5; drug2_c = 7; %Uncomment for ind_mex_dof


female_id = [2002, 2003, 2004, 2005, 2006, 2012, 2013, 2014, 2016];
male_id = [2001, 2007, 2008, 2009, 2010, 2011, 2016, 2017, 2018, 2019, 2020, 2021, 2022];

rel_qtc_male = nan(13, 14); rel_qrs_male = nan(13, 14); rel_tpeakend_male = nan(13, 14); rel_twaveamp_male = nan(13, 14);
male_master_qtc = nan(13, 14); male_master_qrs = nan(13, 14); male_master_tpeakend = nan(13, 14); male_master_twaveamp = nan(13, 14); male_master_drug_conc1 = nan(13,14); male_master_drug_conc2 = nan(13,14); 
female_master_qtc = nan(9, 14); female_master_qrs = nan(9, 14); female_master_tpeakend = nan(9, 14); female_master_twaveamp = nan(9, 14); female_master_drug_conc1 = nan(9,14); female_master_drug_conc2 = nan(9,14);
translated_master_qtc = nan(9, 14); translated_master_qrs = nan(9, 14); translated_master_tpeakend = nan(9, 14); translated_master_twaveamp = nan(9, 14);


for i = 1:13
    % Male data
    male = male_id(i);
    m_ind = find(drug_data.RANDID== male);
    male_id_qtc = table2array(drug_data(m_ind, 21));
    male_id_qrs = table2array(drug_data(m_ind, 16));
    male_id_tpeakend = table2array(drug_data(m_ind, 18));
    male_id_twaveamp = table2array(drug_data(m_ind, 20));
    male_id_drug_conc1_cell = table2array(drug_data(m_ind, drug1_c));
    for kk = 1:length(male_id_drug_conc1_cell); male_id_drug_conc1(kk) = str2double(cell2mat(male_id_drug_conc1_cell(kk))); end

    male_id_drug_conc2_cell = table2array(drug_data(m_ind, drug2_c));
    for kk = 1:length(male_id_drug_conc2_cell); male_id_drug_conc2(kk) = str2double(cell2mat(male_id_drug_conc2_cell(kk))); end


    m_ind_placebo = find(placebo_data.RANDID== male);
    male_id_placebo_qtc = table2array(placebo_data(m_ind_placebo, 21));
    male_id_placebo_qrs = table2array(placebo_data(m_ind_placebo, 16));
    male_id_placebo_tpeakend = table2array(placebo_data(m_ind_placebo, 18));
    male_id_placebo_twaveamp = table2array(placebo_data(m_ind_placebo, 20));

    if (~isempty(male_id_qtc) & ~isempty(male_id_placebo_qtc))

    % Placebo correction
    male_id_qtc(end+1:14)=nan;  male_id_qrs(end+1:14)=nan;   male_id_tpeakend(end+1:14)=nan;    male_id_twaveamp(end+1:14)=nan;
    rel_qtc_male(i,:) = (( (male_id_qtc - male_id_qtc(1)) - (male_id_placebo_qtc - male_id_placebo_qtc(1)) )/male_id_qtc(1) )*100  ;
    rel_qrs_male(i,:) = (( (male_id_qrs - male_id_qrs(1)) - (male_id_placebo_qrs - male_id_placebo_qrs(1)) )/male_id_qrs(1) )*100  ;
    rel_tpeakend_male(i,:) = (( (male_id_tpeakend - male_id_tpeakend(1)) - (male_id_placebo_tpeakend - male_id_placebo_tpeakend(1)) )/male_id_tpeakend(1) )*100 ;
    rel_twaveamp_male(i,:) = (( (male_id_twaveamp - male_id_twaveamp(1)) - (male_id_placebo_twaveamp - male_id_placebo_twaveamp(1)) )/male_id_twaveamp(1) )*100  ;

    qtc_male_abs = (male_id_qtc(1).*(1 + rel_qtc_male(i,:)./100));
    qrs_male_abs = (male_id_qrs(1).*(1 + rel_qrs_male(i,:)./100));
    tpeakend_male_abs = (male_id_tpeakend(1).*(1 + rel_tpeakend_male(i,:)./100));
    twaveamp_male_abs = (male_id_twaveamp(1).*(1 + rel_twaveamp_male(i,:)./100));

    male_master_qtc(i, :) = qtc_male_abs;
    male_master_qrs(i, :) = qrs_male_abs;
    male_master_tpeakend(i, :) = tpeakend_male_abs;
    male_master_twaveamp(i, :) = twaveamp_male_abs;
    male_master_drug_conc1(i,:) = male_id_drug_conc1;
    male_master_drug_conc2(i,:) = male_id_drug_conc2;
   
    end
end
%%
rel_qrs_male_mean = mean(rel_qrs_male, 'omitnan');
rel_qtc_male_mean = mean(rel_qtc_male, 'omitnan');
rel_tpeakend_male_mean = mean(rel_tpeakend_male, 'omitnan');
rel_twaveamp_male_mean = mean(rel_twaveamp_male, 'omitnan');

% Translated_data
Y_simulated_output = zeros(14, 4);
predicted_percentage_change = zeros(14, 4);
for t = 1:14
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

%% Female data
    
for j = 1:9

    female = female_id(j);
    f_ind = find(drug_data.RANDID== female);
    female_id_qtc = table2array(drug_data(f_ind, 21));
    female_id_qrs = table2array(drug_data(f_ind, 16));
    female_id_tpeakend = table2array(drug_data(f_ind, 18));
    female_id_twaveamp = table2array(drug_data(f_ind, 20));
    female_id_drug_conc1_cell = table2array(drug_data(f_ind, drug1_c));
    for kk = 1:length(female_id_drug_conc1_cell); female_id_drug_conc1(kk) = str2double(cell2mat(female_id_drug_conc1_cell(kk))); end
    
    female_id_drug_conc2_cell = table2array(drug_data(f_ind, drug2_c));
    for kk = 1:length(female_id_drug_conc2_cell); female_id_drug_conc2(kk) = str2double(cell2mat(female_id_drug_conc2_cell(kk))); end
     
    f_ind_placebo = find(placebo_data.RANDID== female);
    female_id_placebo_qtc = table2array(placebo_data(f_ind_placebo, 21));
    female_id_placebo_qrs = table2array(placebo_data(f_ind_placebo, 16));
    female_id_placebo_tpeakend = table2array(placebo_data(f_ind_placebo, 18));
    female_id_placebo_twaveamp = table2array(placebo_data(f_ind_placebo, 20));
    
    if ~isempty(female_id_placebo_qtc)
    %Placebo correction
    female_id_qtc(end+1:14)=nan;  female_id_qrs(end+1:14)=nan;   female_id_tpeakend(end+1:14)=nan;    female_id_twaveamp(end+1:14)=nan;
    
    rel_qtc_female = (( (female_id_qtc - female_id_qtc(1)) - (female_id_placebo_qtc - female_id_placebo_qtc(1)) )/female_id_qtc(1) )*100  ;
    rel_qrs_female = (( (female_id_qrs - female_id_qrs(1)) - (female_id_placebo_qrs - female_id_placebo_qrs(1)) )/female_id_qrs(1) )*100  ;
    rel_tpeakend_female = (( (female_id_tpeakend - female_id_tpeakend(1)) - (female_id_placebo_tpeakend - female_id_placebo_tpeakend(1)) )/female_id_tpeakend(1) )*100 ;
    rel_twaveamp_female = (( (female_id_twaveamp - female_id_twaveamp(1)) - (female_id_placebo_twaveamp - female_id_placebo_twaveamp(1)) )/female_id_twaveamp(1) )*100  ;

    qtc_female_abs = (female_id_qtc(1).*(1 + rel_qtc_female./100));
    qrs_female_abs = (female_id_qrs(1).*(1 + rel_qrs_female./100));
    tpeakend_female_abs = (female_id_tpeakend(1).*(1 + rel_tpeakend_female./100));
    twaveamp_female_abs = (female_id_twaveamp(1).*(1 + rel_twaveamp_female./100));

    female_master_qtc(j, :) = qtc_female_abs;
    female_master_qrs(j, :) = qrs_female_abs;
    female_master_tpeakend(j, :) = tpeakend_female_abs;
    female_master_twaveamp(j, :) = twaveamp_female_abs;
    female_master_drug_conc1(j,:) = female_id_drug_conc1;
    female_master_drug_conc2(j,:) = female_id_drug_conc2;
       
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
set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
xticks([0,6,12,18,24 ]);  yticks([80, 90, 100, 110 ]);
%xlabel('Time (h)'); ylabel('QRS dur (ms)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(3); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2016QRS.svg')
% saveas(gcf,'fig_d1/lid_dof_2016QRS.svg')
% saveas(gcf,'fig_d1/mox_dilt_2016QRS.svg')
% saveas(gcf,'fig_d1/mex_dof_2016QRS.svg')

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
set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
xticks([0,6,12,18,24 ]);  yticks([380,420,460,500 ]);
%xlabel('Time (h)'); ylabel('QTc (ms)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(4); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2016QTc.svg')
% saveas(gcf,'fig_d1/lid_dof_2016QTc.svg')
% saveas(gcf,'fig_d1/mox_dilt_2016QTc.svg')
% saveas(gcf,'fig_d1/mex_dof_2016QTc.svg')

figure(5); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(female_master_tpeakend,0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(male_master_tpeakend,0.3,color_m_shade,time_point);
[lineOut3, fillOut3] = stdshade(translated_master_tpeakend,0.3,'k',time_point);
plot(time_point, mean(female_master_tpeakend, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(male_master_tpeakend, 'omitnan'), '.-','Markersize', 12,'linewidth', 1.6, 'Color', color_m_solid)
plot(time_point, mean(translated_master_tpeakend, 'omitnan'), '.-', 'Markersize', 12, 'linewidth', 1.6, 'Color', 'k')
diff_tpeakend_translated = abs((mean(female_master_tpeakend, 'omitnan') - mean(translated_master_tpeakend, 'omitnan'))./mean(female_master_tpeakend, 'omitnan'))';

xlim([-1 25]); ylim([60 150]); 
set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
xticks([0,6,12,18,24 ]);  yticks([60,90,120,150]);
%xlabel('Time (h)'); ylabel('T-peak-end dur. (ms)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(5); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2016Tpeakend.svg')
% saveas(gcf,'fig_d1/lid_dof_2016Tpeakend.svg')
% saveas(gcf,'fig_d1/mox_dilt_2016Tpeakend.svg')
% saveas(gcf,'fig_d1/mex_dof_2016Tpeakend.svg')

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
% saveas(gcf,'fig_d1/dof_2016Twaveamp.svg')
% saveas(gcf,'fig_d1/lid_dof_2016Twaveamp.svg')
% saveas(gcf,'fig_d1/mox_dilt_2016Twaveamp.svg')
% saveas(gcf,'fig_d1/mex_dof_2016Twaveamp.svg')

%{
female_master_drug_conc1(isnan(female_master_drug_conc1)) = 0;   male_master_drug_conc1(isnan(male_master_drug_conc1)) = 0;
figure(7); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(female_master_drug_conc1,0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(male_master_drug_conc1,0.3,color_m_shade,time_point);
hold on 
plot(time_point, mean(female_master_drug_conc1, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(male_master_drug_conc1, 'omitnan'), '.-','Markersize', 12,'linewidth', 1.6, 'Color', color_m_solid)
xlim([-1 25]); 
set(gca, 'xticklabel', {[]});
xticks([0,6,12,18,24 ]);  
%xlabel('Time (h)'); ylabel('T-wave amp. (\muV)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(7); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2016drug.svg')
% saveas(gcf,'fig_d1/lid_dof_2016drug1.svg')
% saveas(gcf,'fig_d1/mox_dilt_2016drug1.svg')
% saveas(gcf,'fig_d1/mex_dof_2016drug1.svg')

female_master_drug_conc2(isnan(female_master_drug_conc2)) = 0;   male_master_drug_conc2(isnan(male_master_drug_conc2)) = 0;
figure(8); hold on; set(gcf, 'color', 'w'); 
[lineOut1, fillOut1] = stdshade(female_master_drug_conc2,0.3,color_f_shade,time_point);
[lineOut2, fillOut2] = stdshade(male_master_drug_conc2,0.3,color_m_shade,time_point);
hold on 
plot(time_point, mean(female_master_drug_conc2, 'omitnan'), '.-', 'Markersize', 12 ,'linewidth', 1.6, 'Color', color_f_solid)
plot(time_point, mean(male_master_drug_conc2, 'omitnan'), '.-','Markersize', 12,'linewidth', 1.6, 'Color', color_m_solid)
xlim([-1 25]); 
set(gca, 'xticklabel', {[]});
xticks([0,6,12,18,24 ]);  
%xlabel('Time (h)'); ylabel('T-wave amp. (\muV)'); 
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.2, 'TickLength', [0.02, 0.02], 'box', 'off', 'tickdir', 'out');
figure(8); set(gcf, 'Units', 'Inches', 'Position', [0 0 4.5 4], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 4])
% saveas(gcf,'fig_d1/dof_2016drug.svg')
% saveas(gcf,'fig_d1/lid_dof_2016drug2.svg')
% saveas(gcf,'fig_d1/mox_dilt_2016drug2.svg')
% saveas(gcf,'fig_d1/mex_dof_2016drug2.svg')

%}
%%

errors_matrix = [diff_qrs_translated diff_qtc_translated diff_tpeakend_translated diff_twaveamp_translated];
errors_matrix(1, :) = [];  %Ignore first point baseline