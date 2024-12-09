% Units
% voltage in mV
% current in uA
% conductance in mS
% resistance in k-oh
% capacitance in uF   I(uF*mV/ms = uA)
% time in ms
% length in um
% concentration in mM

%close all; clear; clc
%Modified by Roshni Feb 2024

%To run the baseline male: runORd_endo_epi(1, 1, 1)
%To run the baseline female: runORd_endo_epi(2, 1, 1)

%To run male population replace par_SA = all_parameters(ii,:); and
%run%runORd_endo_epi(1, 1, 1000)   %Note run smaller populations if computational resources are insufficient
%To run female population replace par_SA = all_parameters(ii,:); and run runORd_endo_epi(2, 1, 1000)

function runORd_endo_epi(gender, start_idx, stop_idx)
tic;

gendertype = gender;  %1: Male; 2: Female

if gendertype ==2
    Ncell = 150 + 40;    Ncell_endo = 95;  Ncell_epi = 95;
 
elseif gendertype ==1
    Ncell = 165 + 40;    Ncell_endo = 102;  Ncell_epi = 103;
end

% cell geometry parameters
L = 100;        % cell length, um
r = 11;         % cell radius, um

%% time parameters
bcl = 1000;      % basic cycle length ms
nbeats = 50;      % number of beats
T = bcl*nbeats;  % Total time

% numerical integration time step 
% use a different time step (dt1) immediately following stimulus
% small dt at stimulus, larger dt during repolarization
dtparams.dt1 = 0.01; %1e-2;      % msec
dtparams.dt1_samp = dtparams.dt1*1;  % sampling time interval
dtparams.dt2 = 0.02;%dtparams.dt1*10;
dtparams.dt2_samp = dtparams.dt2*1;
dtparams.twin = 50; % ms, defines window using dt1 (need to adjust for longer cables and/or slower conduction)

Vthresh = -60;  % activation/repol threshold, mV
%trange = [max(0,(nbeats-4)*bcl-10) nbeats*bcl];  % time range to output values
trange = [0 (nbeats*bcl)];  % time range to output values

%% cell geometry
Aax = 2*pi*r*L; % patch surface area, um^2
Ad = pi*r^2;    % disc surface area, um^2
Atot = 2*Ad + Aax;  % total surface area, um^2

%% bulk extracellular concentrations
Ko = 5.4;                  % mM
Nao = 140;                 % mM
Cao = 1.8;                 % mM

%% model specific parameters

% order is determined by code in fun_name
% INa INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK  IKb INab ICab IpCa
Ncurrents = 14;
scaleI = ones(1, Ncurrents); % ionic current scaling factors
% p is a struct variable consisting of all the scaling factors to create arrays of currents and constants defined in the equation file

p.iina = 1; p.iinal = 2; p.iito = 3; p.iical = 4; p.iikr = 5; p.iiks = 6;
p.iik1 = 7; p.iinaca_i = 8; p.iinaca_ss = 9; p.iinak = 10; p.iikb = 11; p.iinab = 12;
p.iicab = 13; p.iipca = 14;

%% additional scaling factors
p.fSERCA = 1; p.fRyR = 1; p.ftauhL = 1; p.fCaMKa = 1; p. fIleak = 1; p.fJrel = 1;

%% initial conditions
%initial conditions for state variables
x0 = Initial_ORd11;
%X0 is the vector for initial sconditions for state variables

p.Nstate = 41;  % number of state variables, including Vm, per patch

%% stimulus parameters
p.stim_dur = 5;   % ms
p.stim_amp = 50; %5.93*2;    % uA/uF

%% cell typ
p.celltype = [zeros(Ncell_endo,1); ones(Ncell - Ncell_endo,1)];  %endo = 0, epi = 1, M = 2

p.L = L; p.r = r;  % um
% extracellular concentrations
p.K_o = Ko;                  % mM
p.Na_o = Nao;                 % mM
p.Ca_o = Cao;                 % mM
% cleft ionic concentration parameters
p.cleft_conc = [p.Na_o; p.K_o; p.Ca_o];
zvec = [1; 1; 2];  % charge valence

pcable = p; dtparams_cable = dtparams;

%% initial conditions in each cell
% load('Gmat_end_99.mat');
% load('v_end_99.mat');
% pcable.phii_0 = v_end_99; %V after 100s %x0(1);  % Initial V
% g0_vec = Gmat_end_99; %Gmat after 100s
% pcable.g0_vec = g0_vec;

%pcable.phii_0 = x0(1);
pcable.indstim = ismember((1:Ncell)',1:5);  %Stimulate first 5 cells
pcable.ind_tau_ip = 1; pcable.ind_tau_im = 1; pcable.tau_Nai_mat = 1; pcable.tau_Cai_mat = 1; pcable.tau_Ki_mat = 1;

%g0_ vec consists of the initial conditions of all the 50 cells set in one array
g0_vec = zeros(Ncell*(pcable.Nstate-1),1); 
for i = 1:pcable.Nstate-1
    g0_vec(1 + Ncell*(i-1):Ncell*i) = x0(i+1);
end
pcable.g0_vec = g0_vec;
pcable.phii_0 = x0(1);
%% Function generic_cable returs the following outputs
%ts: time vector
%Vm: Voltage of 50 cells across time
%Gmat: State variables across time (other than Vm)
%tup: Upstroke time of 50 cells
%trepol: Repolarization time of 50 cells
%cv_est: Conduction velocity
% p struct variable
% gap junctional coupling

imp_data = importdata('parameter_matrix_1000.mat');
all_parameters = imp_data.all_parameters;

for ii = start_idx: stop_idx

par_SA = ones(1,17); %ones(1,17);%all_parameters(ii,:);

ggap_endo_male = 3400*10^(-6).*par_SA(17);   %mS
if gendertype ==1
    ggap = ones(Ncell-1,1);      % mS 350
    ggap(1: Ncell_endo-1) = ggap_endo_male;
    ggap(Ncell_endo:end) = ggap_endo_male*0.94;
   
elseif gendertype ==2
    ggap = ones(Ncell-1,1);      % mS 350
    ggap(1: Ncell_endo-1) = ggap_endo_male*0.68;
    ggap(Ncell_endo:end) = ggap_endo_male*0.61;
else
    ggap = ones(Ncell-1,1).*ggap_endo_male;
end


APfirst_gKr = 0.04; APlast_gKr = 0.0744; % 

pcable.GKr_male = linspace(APfirst_gKr, APlast_gKr,Ncell )';
pcable.GKr_female =  linspace(0.80*APfirst_gKr, 0.80*APlast_gKr, Ncell )';

[ts, Vm, Gmat, cv_est, pout] =  cable_ORd11_gender(pcable, Ncell, ggap, bcl, scaleI, dtparams_cable, T, Vthresh, trange, gendertype, par_SA); % CV in cm/s

%parsave(Vm, cv_est, gendertype, bcl, ii)

disp(ii)
disp(cv_est)

end
ts_sample = ts(1:5:end,1);
Tname = sprintf('ts_gen%d_Ncell%d_endo%d_epi%d_CL%d.mat', gendertype, Ncell, Ncell_endo, Ncell_epi, bcl); 
%save(Tname, 'ts_sample');

toc;

%% calculate ionic currents
% Iall_mat: all ionic currents matrix (uA), Ncurrents*Ncell x Ns (order as defined in ionic model): grouped by membrane patches, 
% [Current 1 (Cell 1),..., Current Nc (Cell 1),
%  Current 1 (Cell 2),..., Current Nc (Cell 2),..., 
%  Current 1 (Cell Np),..., Current Nc (Cell Ncell)]
% Ncurrents = length(scaleI);
% Sb = [p.cleft_conc(1) p.cleft_conc(2) p.cleft_conc(3)];  % mM, bulk extracellular concentrations
% S_bulk = zeros(3*Ncell, 1);               % cleft ion concentrations
% S_bulk(1:Ncell) = Sb(1); S_bulk(Ncell+1:2*Ncell) = Sb(2); S_bulk(2*Ncell+1:3*Ncell) = Sb(3);
% 
% Iall_mat = nan(Ncurrents*pout.Npatches, length(ts));
% for i = 1:length(ts)
%     [~, ~, ~, ~, Iall_mat(:,i)] = model_ORd11_cable_yang2(ts(i), [Vm(:,i); Gmat(:,i)] , pout, S_bulk,  gendertype);
% end
% INa_all = Iall_mat(p.iina:Ncurrents:end,:);
% IKr_all = Iall_mat(p.iikr:Ncurrents:end,:);
% IKs_all = Iall_mat(p.iiks:Ncurrents:end,:);

%% plots
% 
figure(1); hold on; set(gcf, 'color', 'w'); 
plot(ts, Vm, 'b', 'LineWidth', 1);
xlabel('Time'); ylabel('Voltage(mV)')
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');

% figure(2); hold on; set(gcf, 'color', 'w'); 
% for n=1:Ncell
%     subplot(5,ceil(Ncell/5),n); hold on;
%     plot(ts, Vm(n,:), 'linewidth', 2)
%     title('Cell ' + string(n))
%     xlabel('Time (ms)')
%     ylabel('Vm (mV)')
%     ylim([-90 50])
% end
% set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',15, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');

%%

function parsave(Vm, cv_est, gendertype, bcl, ii)
    Vm_sample = Vm(:, 1:5:end);
    Vname = sprintf('Vm_gen%d_%d_CL%d.mat',gendertype,ii, bcl); save(Vname, 'Vm_sample');
    cv_cable_name = sprintf('cv_cable_gen%d_%d_CL%d.mat', gendertype, ii, bcl); save(cv_cable_name, 'cv_est');
  
end
end