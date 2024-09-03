clear; close all; clc
color = [0 0 0]; % Black
color_m = [0.968627451 0.305882353 0.839215686]; % Pink
color_f = [0.333333333 0.62745098 0.984313725]; % Blue

failed_1x = [39 56 65 71 94 95];
failed_2x = [2 14 25 34  39  45  52  56  65  68 71  72  73  94  95  97 ];
failed_3x = [2 14 25 34  39  42  45  46 52  56 62 65  68 71  72  73  94  95 96 97 98];
failed_4x = [2 8 14 21 24 25 34  39 41  42  45  46 52  56  61  62 65 67 68 71  72  73  94  95 96 97 98];

addpath('../run_simulated_drug_pseudo_ecgs/ecg_features_drug/')

% load('CV_all_1x_male.mat'); [CV_all_1x_male,~] = removerows(CV_all_1x_male,'ind', failed_1x);
load('QRS_amp_1x_male.mat'); [QRS_amp_1x_male,~] = removerows(QRS_amp_1x_male,'ind', failed_1x);
load('QRS_dur_1x_male.mat'); [QRS_dur_1x_male,~] = removerows(QRS_dur_1x_male,'ind', failed_1x);
load('QT_int_1x_male.mat'); [QT_int_1x_male,~] = removerows(QT_int_1x_male,'ind', failed_1x);
load('ST_avg_1x_male.mat'); [ST_avg_1x_male,~] = removerows(ST_avg_1x_male,'ind', failed_1x);
load('T_peakend_dur_1x_male.mat'); [T_peakend_dur_1x_male,~] = removerows(T_peakend_dur_1x_male,'ind', failed_1x);
load('T_wave_amp_1x_male.mat'); [T_wave_amp_1x_male,~] = removerows(T_wave_amp_1x_male,'ind', failed_1x);
load('theta_T_1x_male.mat'); [theta_T_1x_male,~] = removerows(theta_T_1x_male,'ind', failed_1x);
load('Twave_dur_1x_male.mat'); [Twave_dur_1x_male,~] = removerows(Twave_dur_1x_male,'ind', failed_1x);

% load('CV_all_1x_female.mat'); [CV_all_1x_female,~] = removerows(CV_all_1x_female,'ind', failed_1x);
load('QRS_amp_1x_female.mat'); [QRS_amp_1x_female,~] = removerows(QRS_amp_1x_female,'ind', failed_1x);
load('QRS_dur_1x_female.mat'); [QRS_dur_1x_female,~] = removerows(QRS_dur_1x_female,'ind', failed_1x);
load('QT_int_1x_female.mat'); [QT_int_1x_female,~] = removerows(QT_int_1x_female,'ind', failed_1x);
load('ST_avg_1x_female.mat'); [ST_avg_1x_female,~] = removerows(ST_avg_1x_female,'ind', failed_1x);
load('T_peakend_dur_1x_female.mat'); [T_peakend_dur_1x_female,~] = removerows(T_peakend_dur_1x_female,'ind', failed_1x);
load('T_wave_amp_1x_female.mat'); [T_wave_amp_1x_female,~] = removerows(T_wave_amp_1x_female,'ind', failed_1x);
load('theta_T_1x_female.mat'); [theta_T_1x_female,~] = removerows(theta_T_1x_female,'ind', failed_1x);
load('Twave_dur_1x_female.mat'); [Twave_dur_1x_female,~] = removerows(Twave_dur_1x_female,'ind', failed_1x);


% X_test = [CV_all_1x_male, QRS_amp_1x_male, QRS_dur_1x_male, QT_int_1x_male, abs(ST_avg_1x_male), T_peakend_dur_1x_male, T_wave_amp_1x_male, theta_T_1x_male, Twave_dur_1x_male];
% Y_test= [CV_all_1x_female, QRS_amp_1x_female, QRS_dur_1x_female, QT_int_1x_female, abs(ST_avg_1x_female), T_peakend_dur_1x_female, T_wave_amp_1x_female, theta_T_1x_female, Twave_dur_1x_female];
% output_names = {'CV', 'QRSamp', 'QRSdur', 'QTint', 'STavg',...
%     'Tpeakend dur', 'Tawave amp', 'theta Twave', 'Twave dur'};

X_test = [QRS_dur_1x_male, QT_int_1x_male, T_peakend_dur_1x_male, T_wave_amp_1x_male];
Y_test= [QRS_dur_1x_female, QT_int_1x_female, T_peakend_dur_1x_female, T_wave_amp_1x_female];

output_names = {'QRSdur', 'QTint', 'Tpeakend dur', 'Tawave amp'};
test_count = length(QRS_dur_1x_male);

addpath('../../build_regression_model/')
load('regression_Blasso4.mat')
N_outputs_Y = length(output_names);
%% Plotting options

plot_matrix = 1;
plot_fitting = 1;
plot_validation = 1;
flag_plot_application = 1;

male2female = 1;

%%

disp('▣ Validation...')
N_val = test_count;
N_figures = ceil(N_outputs_Y/10);
predicted_outputs = Y_test*0;
for i = 1:N_val
    cell_index = i;
    inputs = X_test(cell_index,:);
    x = log(inputs);
    xz = (x-mean(X_log))./std(X_log);

    yz = xz*Blasso;
    %    yz = xz;
    y = yz.*std(Y_log)+mean(Y_log);
    predicted_outputs(i,:) = exp(y);
end

Y_pred_1x_lasso = predicted_outputs;

% with R^2 = 0.75, the model explains approximately 75% of the variability in the predicted variable
R2ord = zeros(1,N_outputs_Y);
R2adj = zeros(1,N_outputs_Y);
rmse_val = zeros(1,N_outputs_Y);
for i = 1:N_outputs_Y
    mdl = fitlm(Y_test(:,i),predicted_outputs(:,i));
    R2ord(i) = mdl.Rsquared.Ordinary;
    R2adj(i) = mdl.Rsquared.Adjusted;
    rmse_val(i) = mdl.RMSE;
end
R2ord;
R2adj;
R2 = R2adj;
avg_R2_val = mean(R2);

% Residual Standard Deviation
% Standard deviation for a normal distribution, centered on the predicted regression line,
% representing the distribution of actually observed values
oSD_val = zeros(1,N_outputs_Y);
rSD_val = zeros(1,N_outputs_Y);
for i = 1:N_outputs_Y
    oSD_val(i) = std(Y_test(:,i));
    rSD_val(i) = sqrt(sum((predicted_outputs(:,i) - Y_test(:,i) ).^2) / (N_val-2));
end

if plot_validation == 1
    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(3,4,i),hold on
                % Plot data points
                factor = 1;
                %if (output_index == 1) && (i == 6 || i == 7), factor = 1e6; end
                %if (output_index == 2) && (i == 3 || i == 4), factor = 1e6; end
                plot(factor*Y_test(:,dex1),factor*predicted_outputs(:,dex1),'Marker','o','LineStyle','none')%,'Color',color)
                xlabel(['Actual ',output_names{dex1}])
                ylabel(['Predicted ',output_names{dex1}])
                title(['R^2 = ',num2str(R2(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',14)
                % Plot identity line
                ylim_ind = get(gca,'ylim') ;
                xlim_ind = get(gca,'xlim') ;
                minpoint = min([ylim_ind(1),xlim_ind(1)]);
                maxpoint = max([ylim_ind(2),xlim_ind(2)]);
                plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
                xlim([minpoint, maxpoint])
                ylim([minpoint, maxpoint])
                dex1 = dex1+1;
            end
        end
    end

    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(3,4,i),hold on
                % Plot data points
                factor = 1;
                plot(factor*predicted_outputs(:,dex1),(Y_test(:,dex1)-predicted_outputs(:,dex1))/rSD_val(dex1),'Marker','o','LineStyle','none');
                xlabel(['Predicted ', output_names{dex1}])
                ylabel('Standardized residuals')
                title(['rSD/oSD = ',num2str(rSD_val(dex1)/oSD_val(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',14)
                % Plot identity line
                xlim_ind = get(gca,'xlim') ;
                plot([xlim_ind(1), xlim_ind(2)],[0,0],'--k')
                xlim([xlim_ind(1), xlim_ind(2)])

                dex1 = dex1+1;
            end
        end
    end

    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    subplot(2,2,1),bar(R2)%,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names)
    set(gca,'XLim',[0 N_outputs_Y+1])
    ylim([0 1])
    rotateXLabels(gca(), 90)
    title('R^2 values')

    % Add residuals also in this population
    subplot(2,2,2),bar(rSD_val)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD values')

    subplot(2,2,4),bar(oSD_val)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('oSD values')

    subplot(2,2,3),bar(rSD_val./oSD_val)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD/oSD values')

end