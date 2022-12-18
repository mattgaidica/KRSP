startSample = 60*60*53;
useSamples = 60*60*24; % hours
useRange = startSample:startSample+useSamples-1;
odba = T.odba(useRange);
temp = T.temp(useRange);

close all
ff(1200,400);
plot(odba,'k');
yyaxis right;
plot(temp,'r-');

nestSense = temp-smoothdata(normalize(odba),'gaussian',360);
hold on;
plot(nestSense,'b-');