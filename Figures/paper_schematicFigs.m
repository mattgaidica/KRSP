load('/Users/matt/Box Sync/KRSP Axy Data/Temp/BJD5_23867_JUNE-JULY2019_BJD5_1_nest_trim_nest__20190603.mat');

color1 = 'k';
color2 = [254/255,174/255,0];
lw = 2;

close all
ff(1400,800);
subplot(211);
data1 = T.odba(2280:5000);
plot(data1.^2,'color',color2,'linewidth',lw);
ylim([-5 30]);
xlim([-100 numel(data1)+100]);

subplot(212);
data2 = T.odba(15307:18027);
plot(data2,'color',color1,'linewidth',lw);
ylim([-1 6]);
xlim([-100 numel(data2)+100]);

%% night
op = 1;
ff(1400,800);
subplot(211);
data1 = T.odba(2280:5000);
plot(smoothdata((-data1).*rand(size(data1)),'gaussian',30),'color',[color2,op],'linewidth',lw);
ylim([-1.5 0.5]);
xlim([-100 numel(data1)+100]);

subplot(212);
data2 = T.odba(15307:18027);
plot(smoothdata(-data2.*rand(size(data2)),'gaussian',30),'color',[color1,op],'linewidth',lw);
ylim([-1.5 0.5]);
xlim([-100 numel(data2)+100]);