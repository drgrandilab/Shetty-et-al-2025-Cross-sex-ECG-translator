load('drug_data.mat');

%change depending on EPTC to 1x, 2x 3x or 4x


drug_effective = drug_1x;
effect = '1x';
gender = 1;

for drug_no = 1: 99
    IKr_drug_scale = drug_effective(drug_no, 1);
    INaL_drug_scale = drug_effective(drug_no, 2);
    ICaL_drug_scale = drug_effective(drug_no, 3);
    INa_drug_scale = drug_effective(drug_no, 4);
    Ito_drug_scale = drug_effective(drug_no, 5);
    IK1_drug_scale = drug_effective(drug_no, 6);
    IKs_drug_scale = drug_effective(drug_no, 7);

    runORd_endo_epi_drug(gender, drug_no, effect, IKr_drug_scale, INaL_drug_scale, ICaL_drug_scale, INa_drug_scale, Ito_drug_scale, IK1_drug_scale, IKs_drug_scale);

end