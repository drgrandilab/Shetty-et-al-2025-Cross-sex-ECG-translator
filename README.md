# Cross-sex-ECG-translator-Roshni

### Folders
⋅⋅* male_female_cables_pseudo_ecgs: code for running baseline male and female 1D models and to simulate male and female populations. The folder also contains code to obatin the pseudo-ECGs.
⋅⋅* compute_ECG_features: code to calculate pseudo-ECG features. The feature values for the male and female populations have been stored as .mat.
⋅⋅* build_regression_model: lasso regression to build 'bcross', matrix comprising regression coefficients.
⋅⋅* translate_simulated_drugs: code to translate ECG features for simulated drug responses for a dataset of 95 drugs.
⋅⋅* clinical_ECG_translation: code to translate the clinical ECG data from Johanessen *et al.* 2014, 2016.

### Dependencies
All simulations were run using MATLAB 2023b.
OS specifications: Ubuntu 22.04.4 LTS; 64 bit
Processor: Intel® Core™ i7-8700 CPU @ 3.20GHz × 12
