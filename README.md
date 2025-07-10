# Shetty-et-al-2025-Cross-sex-ECG-translator

Shetty R, Morotti S, Sobota V, Bayer JD, Ni H & Grandi E. Development and Clinical Validation of a Cross-Sex Translator of ECG Drug Responses. JACC Clin Electrophysiol; DOI: 10.1016/j.jacep.2025.05.015. [Link](https://www.jacc.org/doi/10.1016/j.jacep.2025.05.015)

### Folders
+ male_female_cables_pseudo_ecgs: code for running baseline male and female 1D models and to simulate male and female populations. The folder also contains code to obatin the pseudo-ECGs.
+ compute_ECG_features: code to calculate pseudo-ECG features. The feature values for the male and female populations have been stored as .mat.
+ build_regression_model: lasso regression to build the 'bcross' matrix comprising regression coefficients.
+ translate_simulated_drugs: code to translate ECG features for simulated drug responses for a dataset of 98 formulations of compounds (including anti-arrhythmic and other miscellaneous pharmacological agents)
+ clinical_ECG_translation: code to translate the clinical ECG data from Johanessen *et al.* 2014, 2016.

### Dependencies
+ All simulations were run using MATLAB 2023b.
+ OS specifications: Ubuntu 22.04.4 LTS; 64 bit
+ Processor: Intel® Core™ i7-8700 CPU @ 3.20GHz × 12

Clinical data used in the study were obtained from Johannesen *et al.* 2014 & 2016 freely available on PhysioNet
1.	Johannesen L, Vicente J, Mason JW, et al. Differentiating Drug-Induced Multichannel Block on the Electrocardiogram: Randomized Study of Dofetilide, Quinidine, Ranolazine, and Verapamil. Clinical Pharmacology & Therapeutics. 2014;96:549–558.

2.	Johannesen L, Vicente J, Mason J, et al. Late sodium current block for drug-induced long QT syndrome: Results from a prospective clinical trial. Clinical Pharmacology & Therapeutics. 2016;99:214–223.

3.	Goldberger AL, Amaral LA, Glass L, et al. PhysioBank, PhysioToolkit, and PhysioNet: components of a new research resource for complex physiologic signals. Circulation. 2000;101:E215-220.
