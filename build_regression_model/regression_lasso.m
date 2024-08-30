clear; close all; clc
addpath('/media/roshni/DATADRIVE01/Roshni/Cross-sex-Translator-RS/compute_ECG_features/')
color = [0 0 0]; % Black
color_m = [0.968627451 0.305882353 0.839215686]; % Pink
color_f = [0.333333333 0.62745098 0.984313725]; % Blue

load('weird_ecg_Jan.mat');  %Abnormal ECGs in population

% load('CV_all_male.mat'); [CV_all_male,~] = removerows(CV_all_male,'ind', weird_ecg_Jan);
load('QRS_amp_male.mat'); [QRS_amp_male,~] = removerows(QRS_amp_male,'ind', weird_ecg_Jan);
load('QRS_dur_male.mat'); [QRS_dur_male,~] = removerows(QRS_dur_male,'ind', weird_ecg_Jan);
load('QT_int_male.mat'); [QT_int_male,~] = removerows(QT_int_male,'ind', weird_ecg_Jan);
load('ST_avg_male.mat'); [ST_avg_male,~] = removerows(ST_avg_male,'ind', weird_ecg_Jan);
load('T_peakend_dur_male.mat'); [T_peakend_dur_male,~] = removerows(T_peakend_dur_male,'ind', weird_ecg_Jan);
load('T_wave_amp_male.mat'); [T_wave_amp_male,~] = removerows(T_wave_amp_male,'ind', weird_ecg_Jan);
load('theta_T_male.mat'); [theta_T_male,~] = removerows(theta_T_male,'ind', weird_ecg_Jan);
load('Twave_dur_male.mat'); [Twave_dur_male,~] = removerows(Twave_dur_male,'ind', weird_ecg_Jan);

% load('CV_all_female.mat'); [CV_all_female,~] = removerows(CV_all_female,'ind', weird_ecg_Jan);
load('QRS_amp_female.mat'); [QRS_amp_female,~] = removerows(QRS_amp_female,'ind', weird_ecg_Jan);
load('QRS_dur_female.mat'); [QRS_dur_female,~] = removerows(QRS_dur_female,'ind', weird_ecg_Jan);
load('QT_int_female.mat'); [QT_int_female,~] = removerows(QT_int_female,'ind', weird_ecg_Jan);
load('ST_avg_female.mat'); [ST_avg_female,~] = removerows(ST_avg_female,'ind', weird_ecg_Jan);
load('T_peakend_dur_female.mat'); [T_peakend_dur_female,~] = removerows(T_peakend_dur_female,'ind', weird_ecg_Jan);
load('T_wave_amp_female.mat'); [T_wave_amp_female,~] = removerows(T_wave_amp_female,'ind', weird_ecg_Jan);
load('theta_T_female.mat'); [theta_T_female,~] = removerows(theta_T_female,'ind', weird_ecg_Jan);
load('Twave_dur_female.mat'); [Twave_dur_female,~] = removerows(Twave_dur_female,'ind', weird_ecg_Jan);

% good_outputs_X = [CV_all_male, QRS_amp_male, QRS_dur_male, QT_int_male, abs(ST_avg_male), T_peakend_dur_male, T_wave_amp_male, theta_T_male, Twave_dur_male];
% good_outputs_Y= [CV_all_female, QRS_amp_female, QRS_dur_female, QT_int_female, abs(ST_avg_female), T_peakend_dur_female, T_wave_amp_female, theta_T_female, Twave_dur_female];
% 
% output_names = {'CV', 'QRS amp', 'QRS dur', 'QT interval', '|Avg ST_{vol}|',...
%     'T-peak-end dur', 'T-wave amp', 'T-wave angle', 'T-wave dur'};

% good_outputs_X = good_outputs_X(:, 9);
% good_outputs_Y = good_outputs_Y(:, 9); 
% output_names = output_names(9);
good_outputs_X = [QRS_dur_male, QT_int_male, T_peakend_dur_male, T_wave_amp_male];
good_outputs_Y= [QRS_dur_female, QT_int_female, T_peakend_dur_female, T_wave_amp_female];

output_names = {'QRSdur', 'QTint', 'Tpeakend dur', 'Tawave amp'};

N_outputs_Y = length(output_names);
corr_all = corrcoef(good_outputs_X);

train_count = 750;  
test_count= 969 - train_count;

X_test = good_outputs_X(end-test_count+1:end,:);
Y_test = good_outputs_Y(end-test_count+1:end,:);

X_train = good_outputs_X(1:end-test_count,:);
Y_train = good_outputs_Y(1:end-test_count,:);

XZ = zscore(log(X_train)); YZ = zscore(log(Y_train));
X_log = log(X_train); Y_log = log(Y_train);

%% Plotting options

plot_matrix = 1;
plot_fitting = 1;
plot_validation = 1;
flag_plot_application = 1;

male2female = 1;

%% Red to white colormap

% Create a heatmap for each of scaling factors


nColors = 64;
whiteColor = [1 1 1];
redColor = [1 0 0];
blueColor = [0 0 1];

% Linearly interpolate between blue and white (reversed)
blueToWhite = [linspace(blueColor(1), whiteColor(1), nColors/2); ...
               linspace(blueColor(2), whiteColor(2), nColors/2); ...
               linspace(blueColor(3), whiteColor(3), nColors/2)]';

% Linearly interpolate between white and red (reversed)
whiteToRed = [linspace(whiteColor(1), redColor(1), nColors/2); ...
              linspace(whiteColor(2), redColor(2), nColors/2); ...
              linspace(whiteColor(3), redColor(3), nColors/2)]';

redWhiteBlueColormap = [blueToWhite; whiteToRed];
colormap(redWhiteBlueColormap);


%% Fitting

% [B4_all, stats] = lasso(X_train, Y_train(:,4), 'CV', 10);
% lassoPlot(B4_all, stats, 'PlotType', 'CV');
% 
% B4_lasso = [stats.Intercept(stats.Index1SE); B4_all(:, stats.Index1SE)];
% 
% Y_pred4 = X_test*B4_all(:, stats.Index1SE) + stats.Intercept(stats.Index1SE);


% [B4_all, stats4] = lasso(X_train_Z, Y_train_Z(:,4), 'CV', 10);
% lassoPlot(B4_all, stats4, 'PlotType', 'CV');
% 
% B4_lasso = [B4_all(:, stats4.Index1SE)];  B4_intercept = stats4.Intercept(stats4.Index1SE);
% 
% Y_pred4 = X_test_Z*B4_lasso;
n_features = length(output_names);
Blasso = nan(n_features, n_features); intercept = nan(n_features,1);
for feature = 1:n_features
    [B_all, stats] = lasso(XZ, YZ(:,feature), 'CV', 10);
    %lassoPlot(B_all, stats, 'PlotType', 'CV');
    B = B_all(:, stats.Index1SE);
    Blasso(:, feature) = B;
    intercept(feature, :) = stats.Intercept(stats.Index1SE);
end

Yz = XZ*Blasso;
Yhat = Yz.*std(Y_log)+mean(Y_log);

% Goodness of fit - R^2
% Fraction of the variance in the dependent variable which is explained by the model
% Calculate agreement of values predicted by regression (Yhat = Bpls*X) with original outputs (Y)
% Assessment on log-transformed values
SSYT = sum((Y_log-ones(train_count,1)*mean(Y_log)).^2);
SSYR = sum((Yhat-ones(train_count,1)*mean(Y_log)).^2);
R2each = SSYR./SSYT;
avg_R2_fit = mean(R2each);

% Assessment on (normal) values
R2ord_fit = zeros(1,N_outputs_Y);
R2adj_fit = zeros(1,N_outputs_Y);
rmse_fit = zeros(1,N_outputs_Y);
for i = 1:N_outputs_Y
    mdl = fitlm(exp(Y_log(:,i)),exp(Yhat(:,i)));
    R2ord_fit(i) = mdl.Rsquared.Ordinary;
    R2adj_fit(i) = mdl.Rsquared.Adjusted;
    rmse_fit(i) = mdl.RMSE;
end
R2ord_fit;
R2adj_fit;

R2_fit = R2adj_fit; % Values plotted in figures
avg_R2calc_fit = mean(R2_fit);

% Residual Standard Deviation
% Standard deviation for a normal distribution, centered on the predicted regression line,
% representing the distribution of actually observed values
oSD = zeros(1,N_outputs_Y);
rSD = zeros(1,N_outputs_Y);
for i = 1:N_outputs_Y
    oSD(i) = std(exp(Y_log(:,i)));
    rSD(i) = sqrt(sum((exp(Yhat(:,i)) - exp(Y_log(:,i)) ).^2) / (train_count-2));
end

% Plot regression coefficients
if plot_matrix == 1
    figure; set(gcf,'color','w')
    imagesc(Blasso'); colormap jet;
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'YDir','normal')
    title('Regression coefficients');
    if male2female == 1
        xlabel('Outputs Male');
        ylabel('Outputs Female');
    else
        xlabel('Outputs Female');
        ylabel('Outputs Male');
    end
    set(gca,'XTick',(1:N_outputs_Y))
    set(gca,'XTickLabel',output_names)
    set(gca,'YTickLabel',output_names)
    set(gca,'YTick',(1:N_outputs_Y))
    rotateXLabels(gca(), 90)
    colorbar

    figure; set(gcf,'color','w')
    h = heatmap(output_names, output_names, corr_all); colormap jet;
    h.Title = 'Correlation coefficients';
    colorbar

    figure(25); set(gcf,'color','w')
    h = heatmap(output_names, output_names, Blasso'); colormap(redWhiteBlueColormap);
    h.CellLabelFormat = '%.2f';
    %h.Title = 'B (lasso)'; 
    h.FontSize = 18; 
    clim([-1, 1])
    axp = struct(h);       %you will get a warning
    axp.Axes.XAxisLocation = 'top';
    colorbar
    set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',24, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
    figure(25); set(gcf, 'Units', 'Inches', 'Position', [0 0 10 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
    %f = gcf; exportgraphics(f, 'Bcross.png', 'Resolution', 300)


end

% Scatter plots
N_figures = ceil(N_outputs_Y/10);

if plot_fitting == 1
    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(3,3,i),hold on
                % Plot data points
                factor = 1;
                plot(factor*exp(Y_log(:,dex1)),factor*exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color)
                xlabel(['Actual ',output_names{dex1}])
                ylabel(['Predicted ',output_names{dex1}])
                title(['R^2 = ',num2str(R2_fit(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',12)
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

    % The error terms are assumed to be:
    % 1) Normally distribuited
    % 2) Homoscedatic (same variance at every X)
    % 3) Independent
    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(3,3,i),hold on
                % Plot data points
                factor = 1;
                plot(factor*exp(Yhat(:,dex1)),(exp(Y_log(:,dex1))-exp(Yhat(:,dex1)))/rSD(dex1),'Marker','o','LineStyle','none','Color',color);
                xlabel(['Predicted ', output_names{dex1}])
                ylabel('Standardized residuals')
                title(['rSD/oSD = ',num2str(rSD(dex1)/oSD(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',12)
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
    subplot(2,2,1),bar(R2_fit,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names)
    set(gca,'XLim',[0 N_outputs_Y+1])
    ylim([0 1])
    rotateXLabels( gca(), 90)
    title('R^2 values')

    subplot(2,2,2),bar(rSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD values')

    subplot(2,2,4),bar(oSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('oSD values')

    subplot(2,2,3),bar(rSD./oSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD/oSD values')
end

%% Validation of the regression model

disp('▣ Validation...')
N_val = test_count;

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
legend_names = output_names;
% {'CV (cm/s)', 'QRS amp. (a.u)', 'QRS dur. (ms)', 'QT interval (ms)', '|Avg ST_{vol}(a.u)|', 'T-peak-end dur. (ms)', 'T-wave amp. (a.u)', 'T-wave angle (\circ)', 'T-wave dur. (ms)'};

if plot_validation == 1
    dex1 = 1;
    for figdex = 1:N_figures
        figure(30)
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:9
            if dex1 <= N_outputs_Y
                h = subplot(3,3,i); hold on
                % Plot data points
                factor = 1;
                %if (output_index == 1) && (i == 6 || i == 7), factor = 1e6; end
                %if (output_index == 2) && (i == 3 || i == 4), factor = 1e6; end
                hold on; plot(nan, nan, '.','color', [1 1 1])
                plot(factor*Y_test(:,dex1),factor*predicted_outputs(:,dex1),'Marker','o','LineStyle','none', 'Color', color_f)%,'Color',color)
                %xlabel(['Actual ',output_names{dex1}])
                %ylabel(['Predicted ',output_names{dex1}])
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
                f=get(gca,'Children');
                legend([f(3)],legend_names{i}); legend boxoff 
            end
        end
    end
    set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
    figure(30); set(gcf, 'Units', 'Inches', 'Position', [0 0 12 10], 'PaperUnits', 'Inches', 'PaperSize', [12, 10])
    %f = gcf; exportgraphics(f, 'validation.png', 'Resolution', 300)

    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(3,3,i),hold on
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

% save regression_Blasso4 Blasso X_train Y_train X_log Y_log output_names


