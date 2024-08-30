#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 12:12:27 2023

@author: roshnishetty
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

ecg_data_master = pd.read_csv('Johannesen_Strauss2014.csv')
column_names_master = ecg_data_master.columns


#%%

ecg_data_c1 = ecg_data_master.drop(columns=['EGREFID', 'AGE', 'HGHT', 'WGHT', 'SYSBP', 'DIABP', 'RACE', 'ETHNIC', 'TPEAKTPEAKP'])

ecg_data_c1['QTc'] = ecg_data_c1['QT']/np.cbrt(ecg_data_c1['RR'])
ecg_data_c1['EXDOSE'].fillna(0, inplace=True)


time_point = [-0.5, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 12, 14, 24]

qtc_all_mean =[];                 qtc_all_std =[];
qt_all_mean =[];                  qt_all_std =[];
qrs_all_mean =[];                 qrs_all_std =[];
blood_conc_all_mean = [];         blood_conc_all_std = []; 

qtc_f_mean =[]; qtc_m_mean =[];                  qtc_f_std =[]; qtc_m_std =[];
qt_f_mean =[]; qt_m_mean =[];                    qt_f_std =[]; qt_m_std =[];
blood_conc_f_mean = []; blood_conc_m_mean = [];  blood_conc_f_std = []; blood_conc_m_std = [];


#drug = 500; drug_name = 'Dofetilide' ;
#drug = 400; drug_name = 'Quinidine' ;
#drug = 1500; drug_name = 'Ranolazine' ;
#drug = 120; drug_name = 'Verapamil' ;
drug = 0; drug_name = 'Placebo' ;

#RR_all = ecg_data_c1.query('  EXDOSE == @drug ')['RR']
#qt_all = ecg_data_c1.query(' EXDOSE == @drug ')['QT']
RR_all = ecg_data_c1['RR']
qt_all = ecg_data_c1['QT']
mask = ~np.isnan(RR_all) & ~np.isnan(qt_all)
slope, intercept, r_value, p_value, std_err = stats.linregress(RR_all[mask], qt_all[mask])

#QTc= QT + slope of linear regression (1-RR)
ecg_data_c1['QTc_new'] = ecg_data_c1['QT'] + (slope*1000)*(1- RR_all/1000)

#%%
for k in range(0,16):

 
    qtc_all = ecg_data_c1.query(' TPT == @time_point[@k] & EXDOSE == @drug ')['QTc_new']
    blood_conc_all = ecg_data_c1.query(' TPT == @time_point[@k] & EXDOSE == @drug  ')['PCSTRESN']
    qt_all = ecg_data_c1.query(' TPT == @time_point[@k] & EXDOSE == @drug  ')['QT']
  
    qtc_f = ecg_data_c1.query(' SEX == "F" & TPT == @time_point[@k] & EXDOSE == @drug  ')['QTc_new']
    blood_conc_f = ecg_data_c1.query(' SEX == "F" & TPT == @time_point[@k] & EXDOSE == @drug  ')['PCSTRESN']
    qt_f = ecg_data_c1.query(' SEX == "F" & TPT == @time_point[@k] & EXDOSE == @drug  ')['QT']
  
    qtc_m = ecg_data_c1.query(' SEX == "M" & TPT == @time_point[@k] & EXDOSE == @drug  ')['QTc_new']
    blood_conc_m = ecg_data_c1.query(' SEX == "M" & TPT == @time_point[@k] & EXDOSE == @drug  ')['PCSTRESN']  
    qt_m = ecg_data_c1.query(' SEX == "M" & TPT == @time_point[@k] & EXDOSE == @drug  ')['QT']


    qtc_f_mean.append(np.mean(qtc_f));                    qtc_m_mean.append(np.mean(qtc_m))
    qt_f_mean.append(np.mean(qt_f));                      qt_m_mean.append(np.mean(qt_m))
    blood_conc_f_mean.append(np.mean(blood_conc_f));      blood_conc_m_mean.append(np.mean(blood_conc_m))    

    qtc_f_std.append(np.std(qtc_f));                      qtc_m_std.append(np.std(qtc_m))
    qt_f_std.append(np.std(qt_f));                        qt_m_std.append(np.std(qt_m))
    blood_conc_f_std.append(np.std(blood_conc_f));        blood_conc_m_std.append(np.std(blood_conc_m))   
    
    qtc_all_mean.append(np.mean(qtc_all));  
    qt_all_mean.append(np.mean(qt_all));  
    blood_conc_all_mean.append(np.mean(blood_conc_all));  
    
    qtc_all_std.append(np.std(qtc_all));  
    qt_all_std.append(np.std(qt_all)); 
    blood_conc_all_std.append(np.std(blood_conc_all));      

#%%

'''
plt.figure(dpi = 500);
plt.figsize=(15, 6)
plt.errorbar(time_point, blood_conc_f_mean, blood_conc_f_std, linewidth=1, linestyle='-', marker='.', color ='b')
plt.errorbar(time_point, blood_conc_m_mean, blood_conc_m_std, linewidth=1, linestyle='-', marker='.', color ='m')
plt.title('Dofetilide')
plt.xlabel('Time point (hours)')
plt.ylabel('Blood conc. (pg/ml)')


plt.figure(dpi = 500);
plt.figsize=(15, 6)
plt.errorbar(time_point, blood_conc_f_mean, blood_conc_f_std, linewidth=1, linestyle='-', marker='.', color ='b')
plt.errorbar(time_point, blood_conc_m_mean, blood_conc_m_std, linewidth=1, linestyle='-', marker='.', color ='m')
plt.errorbar(time_point, blood_conc_all_mean, blood_conc_all_std, linewidth=1, linestyle='-', marker='.', color ='k')
plt.title('Dofetilide')
plt.xlabel('Time point (hours)')
plt.ylabel('Blood conc. (pg/ml)')

plt.show()

#%%


plt.figure(dpi = 500);
plt.figsize=(15, 6)
plt.errorbar(time_point, qtc_f_mean, qtc_f_std, linewidth=1, linestyle='-', marker='.', color ='b')
plt.errorbar(time_point, qtc_m_mean, qtc_m_std, linewidth=1, linestyle='-', marker='.', color ='m')
plt.title('Dofetilide')
plt.xlabel('Time point (hours)')
plt.ylabel('QTc. (ms)')


plt.figure(dpi = 500);
plt.figsize=(15, 6)
plt.errorbar(time_point, qtc_f_mean, qtc_f_std, linewidth=1, linestyle='-', marker='.', color ='b')
plt.errorbar(time_point, qtc_m_mean, qtc_m_std, linewidth=1, linestyle='-', marker='.', color ='m')
plt.errorbar(time_point, qtc_all_mean, qtc_all_std, linewidth=1, linestyle='-', marker='.', color ='k')
plt.title('Dofetilide')
plt.xlabel('Time point (hours)')
plt.ylabel('QTc. (ms)')

plt.show()
'''
#%%


plt.figure(dpi = 500);
plt.figsize=(15, 6)
plt.errorbar(time_point, blood_conc_f_mean, blood_conc_f_std, linewidth=1, linestyle='-', marker='.', color ='b')
plt.errorbar(time_point, blood_conc_m_mean, blood_conc_m_std, linewidth=1, linestyle='-', marker='.', color ='m')
plt.errorbar(time_point, blood_conc_all_mean, blood_conc_all_std, linewidth=1, linestyle='-', marker='.', color ='k')
plt.title(drug_name)
plt.xlabel('Time point (hours)')
plt.ylabel('Blood conc. (ng/ml)')


plt.figure(dpi = 500);
plt.figsize=(15, 6)
plt.errorbar(time_point, qtc_f_mean, qtc_f_std, linewidth=1, linestyle='-', marker='.', color ='b')
plt.errorbar(time_point, qtc_m_mean, qtc_m_std, linewidth=1, linestyle='-', marker='.', color ='m')
plt.title(drug_name)
plt.xlabel('Time point (hours)')
plt.ylabel('QTc manual (ms)')

plt.figure(dpi = 500);
plt.figsize=(15, 6)
plt.errorbar(time_point, qt_f_mean, qt_f_std, linewidth=1, linestyle='-', marker='.', color ='b')
plt.errorbar(time_point, qt_m_mean, qt_m_std, linewidth=1, linestyle='-', marker='.', color ='m')
plt.title(drug_name)
plt.xlabel('Time point (hours)')
plt.ylabel('QT. (ms)')

#%%

