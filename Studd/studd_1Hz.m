% dataFile = '/Users/matt/Documents/Data/KRSP/J1_Sep1_2014.csv';
% dataFile = '/Users/matt/Documents/Data/KRSP/II_Sept28_2014_1.csv'; % no
% headers
dataFile = '/Users/matt/Documents/Data/KRSP/K6_Sep1_2014.csv';
% dataFile = '/Users/matt/Documents/Data/KRSP/H2_Sep1_2014.csv';
% dataFile = '/Users/matt/Documents/Data/KRSP/O1_Sep1_2014.csv';
% A = readtable(dataFile);
close all
dates = unique(A.date);

ff(1200,600);
% plot(A.X(n));
% hold on;
% plot(A.Y(n));
% plot(A.Z(n));

% yyaxis right;
% smoothtemp = eegfilt(A.temp(n)',1,0,1/100);
for iDay = 1:numel(dates)
    n = (A.date == dates(iDay));
    subplot(211);
    plot(smoothTemp(A.temp(n),[22.01,22.5]));
    hold on;
    title([num2str(numel(dates)),' days temp']);
    ylabel('temp (C)');
    xlabel('time (s)');
    xlim([1 numel(A.temp(n))]);
    
    subplot(212);
    plot(smooth(A.odba(n),200));
    hold on;
    title([num2str(numel(dates)),' days odba']);
    ylabel('odba (g)');
    xlabel('time (s)');
    xlim([1 numel(A.temp(n))]);
end