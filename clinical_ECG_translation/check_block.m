clear; close all; clc

color_list = {[1, 0.6, 0.6], [0.6, 1, 0.6], [0.6, 0.8, 1], [1, 0.8, 0.6], [0.8, 0.6, 1], [0.6, 1, 1], [1, 0.6, 1]}; 
% Order: IKr, INaL, ICaL, INa, Ito, IK1, IKs
%Dofetilide
% IC50 = [4.9, 753160.4, 260.3, 380.5, 18.8, 394.3 ];
% nh = [0.9, 0.3, 1.2, 0.9, 0.8, 0.8];
% MW = 441.567;%  %g/mol

%Quinidine
% IC50 = [992, 9417, 51592.3, 12329, 3487.4, 39589919, 4898.9 ];
% nh = [0.8, 1.3, 0.6, 1.5, 1.3, 0.4, 1.4];
% MW = 782.94;   %g/mol

%Verapamil HCl
IC50 = [288, 7028, 201.8, nan, 13429.2, 349000000, nan];
nh = [1, 1, 1.1, 1, 0.8, 0.3, 1];
MW = 491.06;   %g/mol

% Ranozaline
% IC50 = [8270, 7884.5, 0, 68774, 0, 0, 36155020];
% nh = [0.9, 0.9, 1, 1.4, 1, 1, 0.5];
% MW = 427.537;   %g/mol

% Diltiazem
% IC50 = [13150, 21868.5, 112.1, 110859, 2820000000];
% nh = [0.9, 0.7, 0.7, 0.7, 0.2];
% MW = 414.519;   %g/mol


%%
k_male = zeros(length(nh),1);
k_female = zeros(length(nh),1);

for i = 1: length(nh)
%nM
nh_i = nh(i);%0.8;%0.9;
IC50_i = IC50(i);
d = [0: 0.0001: 100].*IC50_i;
k = 1./ (1+ (d./IC50_i).^nh_i);

male_conc = 84.0636*1000/MW;  %1425*1000/MW; 
female_conc = 164.3727*1000/MW; %2065*1000/MW;
[min_val_male, idx_male_k] = min(abs(d -male_conc));
[min_val_female, idx_female_k] = min(abs(d -female_conc));

k_male(i) = k(idx_male_k);   k_female(i) = k(idx_female_k); 

figure(1);hold on; set(gcf, 'color', 'w'); 
plot(d, k, '.-', 'LineWidth',2, 'Color', color_list{i});
xlabel('Drug (nM)'); ylabel('k');
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
xlim([0, 350]); %xlim([0, 2700]); 
ylim([0, 1])

% figure(2);hold on; set(gcf, 'color', 'w'); 
% plot(d*MW/1000, k, '.-', 'LineWidth',2); 
% xlabel('Drug (ng/mL)'); ylabel('k');
% set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
end
                     
figure(1); hold on;
xline(male_conc, '-', 'Male', 'Color', 'k' , 'Linewidth', 1.5); 
xline(female_conc, '-', 'Female', 'Color','k', 'Linewidth', 1.5 ); 
legend('IKr', 'INaL', 'ICaL', 'INa', 'Ito', 'IK1', 'IKs') ;

