load('algorithmApp_session.mat');
% close all
ff(1200,400);
plot(sense.tempZ(useRange),'r-');
hold on;
plot(binNestSense(useRange),'k-');

yyaxis right;
plot(normalize(gradient(sense.tempZ(useRange))),'r:');
hold on;
plot(normalize(smoothdata(diff(sense.tempZ(useRange),2),'gaussian',60)),'b:');


%%
[BF, BC] = bimodalitycoeff(alg_temp)
close all
ff(400,300);
histogram(alg_temp);

[BF, BC] = bimodalitycoeff(smOdba)
ff(400,300);
histogram(smOdba);

