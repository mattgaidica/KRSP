if do
    sqkey = readtable('sqkey.txt');
    Tss = makeTss(2014:2020);
    iSq = 'G10_Sep1_2014__20140909.mat';
    T = loadTStruct(iSq,sqkey,Tss);
    do = false;
end

color1 = [0 0 0];
color2 = [254/255,174/255,0];
lw = 2;

close all
ff(1400,800);
subplot(211);
data1 = T.odba((2280:5000)+900);
plot(data1.^2,'color',color2,'linewidth',lw);
ylim([-5 30]);
xlim([-100 numel(data1)+100]);

subplot(212);
data2 = T.odba((15307:18027)+900);
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

%% layering
dR = (15307:18027)+900;
dF = 10;

lw = 5;
colors = magma(10);
% close all; clc
ff(700,600);
plot(dec_data2,'color',color1,'linewidth',lw);
ylim([-1 6]);
xlim([-10 numel(dec_data2)+10]);
hold on;

filtData = smoothdata(T.odba,'loess',numel(T.odba)/(dF));
plot(normalize(decimate(filtData(dR),dF),'range')*2+1,'linewidth',lw,'color',colors(2,:));

filtData = smoothdata(T.odba,'loess',numel(T.odba)/(4*dF));
plot(normalize(decimate(filtData(dR),dF),'range')*2+2,'linewidth',lw,'color',colors(5,:));

filtData = smoothdata(T.odba,'loess',numel(T.odba)/(12*dF));
plot(normalize(decimate(filtData(dR),dF),'range')*2+3,'linewidth',lw,'color',colors(10,:));

%% homeograph
% run with doPlot = 1; >> T = loadTStruct(iSq,sqkey,Tss);
% select top subplot
xlim([min(dR) max(dR)]); % see above
% select left axis
ylim([0 60]);
% select right axis
ylim([-10 10]);
yticks(0);
ylabel('Homeograph (H)');
xticks([]);
xlabel('Time');
set(gca,'fontsize',20)
% zoom in
xlim(17698,17698+100);
