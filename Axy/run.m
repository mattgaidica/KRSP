[BF, BC] = bimodalitycoeff(alg_temp)
close all
ff(400,300);
histogram(alg_temp);

[BF, BC] = bimodalitycoeff(smOdba)
ff(400,300);
histogram(smOdba);

