function [QRS_dur, QRS_amp, QT_int, ST_avg, Twave_dur, T_peakend_dur, Twave_amp, theta_T] = get_ECG_features(time_i, ecg_i, gendertype, num)

try 
    time = time_i./1000 - 45;
    int_index = (time > 3.98 & time < 4.98);
    time_int = time(int_index);
    
    ecg = ecg_i;
    ecg_int = ecg(int_index);
    
    dVm_array = (ecg_int(2:end)-ecg_int(1:end-1))./(time_int(2:end).*1000-time_int(1:end-1).*1000); 
    dVm = [dVm_array(1); dVm_array];
    
    Q_ind = find(abs(time_int - 4) < 0.0001);   Q_point = time_int(Q_ind(2));
    [R_point, R_ind] = max(ecg_int);  
    
    T_time_ind = time_int > 4.1;
    T_time_int = time_int(T_time_ind);
    [T_peak_value, T_peak_ind] = max(ecg_int(T_time_ind));
    
    T_wave_on = T_time_int(T_peak_ind) - 0.09;
    T_wave_off = T_time_int(T_peak_ind) + 0.09;
    
    T_interest_ind = (time_int > T_wave_on & time_int < T_wave_off);
    T_interest = time_int(T_interest_ind);
    
    T_interest_ind_up = (time_int > T_wave_on & time_int < T_time_int(T_peak_ind));
    T_interest_up = time_int(T_interest_ind_up);
    ecg_interest_up = ecg_int(T_interest_ind_up);
    
    T_interest_ind_down = (time_int > T_time_int(T_peak_ind) & time_int < T_wave_off);
    T_interest_down = time_int(T_interest_ind_down);
    ecg_interest_down = ecg_int(T_interest_ind_down);
    
    [slope_max_up, slope_max_up_ind] = max(dVm(T_interest_ind_up));
    [slope_max_down, slope_max_down_ind] = min(dVm(T_interest_ind_down));

    theta_T = atand(slope_max_up*1000);
    
    x_tangent_up = T_interest_up(slope_max_up_ind);
    y_tangent_up = ecg_interest_up(slope_max_up_ind);
    ytanget_line_up = slope_max_up*1000*(time_int - x_tangent_up) + y_tangent_up;
    
    [val0up,idx0up]=min(abs(ytanget_line_up));
    minVal_up= ytanget_line_up(idx0up);
    t_on_time = time_int(idx0up);
    
    x_tangent_down = T_interest_down(slope_max_down_ind);
    y_tangent_down = ecg_interest_down(slope_max_down_ind);
    ytanget_line_down = slope_max_down*1000*(time_int - x_tangent_down) + y_tangent_down;
    
    [val0down,idx0down]=min(abs(ytanget_line_down));
    minVal_down= ytanget_line_down(idx0down);
    t_off_time = time_int(idx0down);
    
    S_time_int =  time_int(R_ind: idx0up);
    S_ecg_int  = ecg_int(R_ind: idx0up);
    amp_S_ind = find(S_ecg_int < 0.01);
    S_point = S_time_int(amp_S_ind(1));
    
    [S_amp,S_idx]=min(abs(time_int - S_point));
    ST_seg = ecg_int(S_idx: idx0up);
    
    
    h1 = figure(); hold on
    yyaxis left
    plot(time_int, ecg_int, 'k'); hold on;
    
    plot(time_int, ytanget_line_up, 'r'); hold on
    plot(time_int, ytanget_line_down, 'r'); hold on
    scatter(Q_point, ecg_int(Q_ind), 'b', 'filled');
    scatter(time_int(R_ind), ecg_int(R_ind), 'r', 'filled');
    scatter(S_point, S_ecg_int(amp_S_ind(1)), 'g', 'filled');
    %scatter(time_int(S_idx), ecg_int(S_idx), 'k', 'filled')
    scatter(t_on_time, minVal_up, 'y', 'filled');
    scatter(T_time_int(T_peak_ind), T_peak_value, 'm', 'filled');
    scatter(t_off_time, minVal_down, 'c', 'filled');
    scatter(T_interest_up(slope_max_up_ind), ecg_interest_up(slope_max_up_ind), '*m');
    scatter(T_interest_down(slope_max_down_ind), ecg_interest_down(slope_max_down_ind), '*m');
    ylim([-0.15 1.3])
    
    yyaxis right
    plot(time_int, dVm, 'b'); hold on
    xlim([3.9 4.7])
    
    title(['Sex:' , num2str(gendertype) , '  Index' , num2str(num)])
    % pause(1)
    % figname = sprintf('/figures_ECG/FIG_%d_%d.png',gendertype, num);
    % saveas(h1,[pwd figname]);
        
    close all
    
    QRS_dur = (S_point - Q_point)*1000;
    QT_int = (t_off_time - Q_point)*1000;
    %ST_midpoint = (minVal_up -  S_ecg_int(amp_S_ind(1)))/2;
    Twave_dur = (t_off_time - t_on_time)*1000;
    QRS_amp =  ecg_int(R_ind);
    Twave_amp= T_peak_value;
    ST_avg = mean(ST_seg); 
    T_peakend_dur = (t_off_time - T_time_int(T_peak_ind))*1000;


catch
    disp(num)
end


end



%%

