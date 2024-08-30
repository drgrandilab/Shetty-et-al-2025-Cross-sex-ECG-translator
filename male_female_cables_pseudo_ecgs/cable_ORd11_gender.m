function [time_final5, phi_final5, Gmat_final5, cv_est, p] =  cable_ORd11_gender(p, Ncell, ggap, bcl, scaleI, dtparams, T, Vthresh, trange, gendertype, par_SA)

% Modified by Roshni Shetty July 2023;
% Code provided by Weinberg Lab, Ohio State University
% https://bme.osu.edu/weinberg-lab-computational-physiology

F_E = 0;
% Units
%
% voltage in mV
% current in uA
% conductance in mS
% resistance in k-ohm
% capacitance in uF    (uF*mV/msec = uA)
% time in msec
% length in um
% concentration in mM

if isempty(trange)
    trange = [0 T];
end

%% time step parameters

% get_time_vector generates a cleverly defined vector - small dt immediately after stimulus, larger dt during repolarization
dt1 = dtparams.dt1; dt2 = dtparams.dt2; 
dt1_samp = dtparams.dt1_samp; dt2_samp = dtparams.dt2_samp;
twin = dtparams.twin;
dt2_new = 0.06;

%ts = get_time_vector(trange, dt1_samp, dt2_samp, twin, bcl);

ts_save = get_time_vector(trange, 0.01, 0.06, twin, bcl);

p.bcl = bcl;       % CL
L = p.L; r = p.r;   % cell dimensions, um

Cm = 1*1e-8;      % membrane capacitance, uF/um^2
Aax = 2*pi*r*L; % patch surface area, um^2
Ad = pi*r^2;    % disc surface area, um^2
Atot = 2*Ad + Aax;  % total surface area, um^2
Ctot = Cm*Atot; % total cell capacitance, uF

p.Ctot = Atot*Cm;   % total cell capacitance, uF

Npatches = Ncell; % number of membrane patches
p.N = Ncell;
p.Npatches = Npatches;
% mS, myoplasmic conductance, scalar


% indexing for phi terms (in mV):
% phi_m,i:               i = 1:N

% indexing for I terms (in uA) and V terms (in mV): (for membrane patches)
% V/I_m,i:            i = 1:N


f_I = scaleI'*ones(1, Npatches); % scaling factors for ionic currents
p.f_I = f_I;

G = zeros(Npatches*(p.Nstate-1),1);  % gating variables
S = zeros(3*Npatches, 1);            % ion concentrations
phi = zeros(Npatches, 1);            % voltages

% initial conditions
phi(:) = p.phii_0;  % potential
Sb = [p.cleft_conc(1) p.cleft_conc(2) p.cleft_conc(3)];  % mM, bulk extracellular concentrations
G(:) = p.g0_vec; % gating variables
S(1:Ncell) = Sb(1);  % extracellular Na+
S(Ncell+1:2*Ncell) = Sb(2);  % extracellular K+
S(2*Ncell+1:3*Ncell) = Sb(3);  % extracellular Ca2+
S_bulk = S;

%phi_mat = nan(length(phi),length(ts));   % Voltage of cells across time
%Gmat = nan(length(G), length(ts));       % Remaining state variables of cells across time

%Save last 5 beats 
ind_last = find(ts_save > (T - 5000));
time_final5 = nan(length(ind_last), 1);%ts_save(ind_last);
phi_final5 = nan(length(phi),length(ind_last));    % Voltage of cells across time
% time_final5 = nan(104167, 1);%ts_save(ind_last);
% phi_final5 = nan(length(phi),104167);    % Voltage of cells across time
% 

Gmat_final5 = 0;%nan(length(G), length(ind_last));   


if trange(1)==0  %At t0
    %phi_mat(:,1) = phi;
    %Gmat(:,1) = G;
    count = 1;
    count2 = 1;
else
    count = 0;
end

tic; %start timer

Vm = phi;

Vm_old = Vm;
beat_num = ones(Npatches,1);
tup = nan(Npatches,1);  %upstroke times
trepol = nan(Npatches,1); %repolarization times


ind_i = 2:Ncell-1;   % indices of axial membrane patch
ind_p = 3:Ncell;
ind_m = 1:Ncell-2;

ti = 0;

while ti < T
    %set dt
    if mod(ti,bcl)<twin
        dt = dt1; dt_samp = dt1_samp;
    else
        dt = dt2; dt_samp = dt2_samp;
    end

    if mod(ti,bcl)<twin
        dt_new = dt1; 
    else
        dt_new = dt2_new ;
    end
    
    p.dt = dt;
    
    
    % update ionic currents, gating variables
    if F_E == 1
        [G_new, Iion, ~, ~, ~] = model_ORd11_cable_alex_ORd2021(ti,[Vm;G],p, S_bulk, gendertype, par_SA);  % 1Male, 2Female
        
        %cable equations
        dV = nan(size(Vm));
        dV(1) = (1/Ctot)*(ggap(1)*(Vm(2)-Vm(1)) - Iion(1));
        dV(2:end-1) = (1/Ctot)*(ggap(1:end-1).*(Vm(ind_m)-Vm(ind_i)) + ggap(2:end).*(Vm(ind_p)-Vm(ind_i)) - Iion(2:end-1));
        dV(end) = (1/Ctot)*(ggap(end)*(Vm(end-1)-Vm(end)) - Iion(end));
        
        phi_new = phi + dt*dV;
        
        %phi_new_dt2 = phi + dt_new*dV;
        
        if any(isnan(Iion))
            disp('Error'); break;
        end
    else
        %Step1
        dV = nan(size(Vm));
        dV(1) = (1/Ctot)*(ggap(1)*(Vm(2)-Vm(1)));
        dV(2:end-1) = (1/Ctot)*(ggap(1:end-1).*(Vm(ind_m)-Vm(ind_i)) + ggap(2:end).*(Vm(ind_p)-Vm(ind_i)));
        dV(end) = (1/Ctot)*(ggap(end)*(Vm(end-1)-Vm(end)));

        phi_new = phi + (dt/2)*dV;

        %Step2
        [G_new, Iion, ~, ~, ~] = model_ORd11_cable_alex_ORd2021(ti,[phi_new;G],p, S_bulk, gendertype, par_SA);  % 1Male, 2Female
        
        if any(isnan(Iion))
           disp('Error Nan Iion'); break;
        end


        phi_new2 = phi_new - dt*Iion/Ctot;
 

        %Step3
        dV(1) = (1/Ctot)*(ggap(1)*(phi_new2(2)-phi_new2(1)));
        dV(2:end-1) = (1/Ctot)*(ggap(1:end-1).*(phi_new2(ind_m)-phi_new2(ind_i)) + ggap(2:end).*(phi_new2(ind_p)-phi_new2(ind_i)));
        dV(end) = (1/Ctot)*(ggap(end)*(phi_new2(end-1)-phi_new2(end)));

        phi_new = phi_new2  + (dt/2)*dV;
        
    end

    Vm = phi_new;
    
    if any(isnan(phi_new))
        disp('Error'); break;
    end
    
    % calculate activation/repolarization times
    [tup, trepol, beat_num] = update_tup_repol(ti, dt, Vm, Vm_old, Vthresh, tup, trepol, beat_num);
    
    Vm_old = Vm;
    
    ti = round(ti + dt, 4);
    %ti_new = round(ti_new + dt_new, 4);
    
    if ti>=trange(1) && ti<=trange(2) && ~mod(ti, dt_samp)
        count = count + 1;
        %phi_mat(:,count) = phi_new;
        %Gmat(:,count) = G_new;
    end
    
    %if ti_new > (trange(2) - 5000) && ~mod(ti_new, dt_new) %ismember(ti, time_final5)
    if ti > (trange(2) - 5000) && ~mod(ti, dt_new) %ismember(ti, time_final5)
        
        phi_final5(:,count2) = phi_new;
        time_final5(count2, 1) = ti;
        %Gmat_final5(:,count2) = G_new;
        count2 = count2 + 1;
    end
    
    
    phi = phi_new;
    G = G_new;
    
end
toc; %Stop timer

% simple conduction velocity calculation
% use middle 50% of cable for CV
i1 = round(.25*Ncell); i2 = round(.75*Ncell);
cv_est = 100*(i2-i1)*(p.L/1000)./(tup(i2,:)-tup(i1,:));  % mm/ms = m/s, 100*m/s = cm/s

% cv_endo = 100*(120-40)*(p.L/1000)./(tup(120,:)-tup(40,:));  % mm/ms = m/s, 100*m/s = cm/s
% cv_epi = 100*(310-210)*(p.L/1000)./(tup(310,:)-tup(210,:));  % mm/ms = m/s, 100*m/s = cm/s
% cv_transmural = 100*(165-155)*(p.L/1000)./(tup(165,:)-tup(155,:));  % mm/ms = m/s, 100*m/s = cm/s
% fprintf('Endo CV: %.2f \n', round(cv_endo(end), 2));
% fprintf('Epi CV:%.2f \n', round(cv_epi(end), 2));
% fprintf('Transmural CV:%.2f \n', round(cv_transmural(end), 2));

end

%% calculate activation/repolarization times
function [tup, trepol, beat_num] = update_tup_repol(ti, dt, Vm, Vm_old, Vthresh, ...
    tup, trepol, beat_num)
iact = find((Vm>Vthresh).*(Vm_old<Vthresh));
for j = 1:length(iact)
    y1 = Vm_old(iact(j)); y2 = Vm(iact(j));
    m = (y2-y1)/dt;
    tup(iact(j),beat_num(iact(j))) = ti - (y1-Vthresh)./m;
    %         save(times_name,'tup','trepol');
end
irepol = find((Vm<Vthresh).*(Vm_old>Vthresh));
for j = 1:length(irepol)
    if ti-tup(irepol(j),beat_num(irepol(j)))>.1
        y1 = Vm_old(irepol(j)); y2 = Vm(irepol(j));
        m = (y2-y1)/dt;
        trepol(irepol(j),beat_num(irepol(j))) = ti - (y1-Vthresh)./m;
        beat_num(irepol(j)) = beat_num(irepol(j)) + 1;
    end
    %         save(times_name,'tup','trepol');
end

end

%% obtain spaced out time vector
function ts = get_time_vector(trange, dt1, dt2, twin, bcl)

ti = trange(1);  count = 1;
ts = nan(1,1e6);
ts(1) = trange(1);
while ti < trange(2)
    if mod(ti,bcl)<twin
        dt = dt1;
    else
        dt = dt2;
    end

    ti = round(ti + dt,4);
    ts(count+1) = ti;
    count = count + 1;
end
ts = ts(1:count);
end

%%


