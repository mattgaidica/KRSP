if do
    T = readtable('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/Emily Studds Axy/April 2015 - Week 1/B10_Apr23_2015.csv');
end

startSample = 60*60*5;
useSamples = 60*60*2; % 1 hour
x = T.X(startSample:startSample+useSamples-1);
y = T.Y(startSample:startSample+useSamples-1);
z = T.Z(startSample:startSample+useSamples-1);
nest = strcmp(T.Nest(startSample:startSample+useSamples-1),'Nest');

colors = lines(5);
t = linspace(0,useSamples/60,numel(x)); % minutes
close all
ff(1200,400);
plot(t,x);
hold on;
plot(t,y);
plot(t,z);
ylabel('raw axy');
yyaxis right;
plot(t,nest,'k');
set(gca,'ycolor','k');
ylim([-1 2]);
yticks([0 1]);
yticklabels({'Out of Nest','In Nest'});
xlim([min(t),max(t)]);
xlabel('Time (min)');

%% dead reckon
d_x = cumsum(diff(x,2));
d_y = cumsum(diff(y,2));
d_z = cumsum(diff(z,2));
t_d = linspace(0,useSamples/60,numel(d_x)); % minutes
close all;
ff(1000,400);
plot(t_d,d_x);
hold on;
plot(t_d,d_y);
plot(t_d,d_z);
ylabel('dead axy');
yyaxis right;
plot(t_d,nest(1:numel(d_x)),'k');
set(gca,'ycolor','k');
ylim([-1 2]);
yticks([0 1]);
yticklabels({'Out of Nest','In Nest'});
xlim([min(t_d),max(t_d)]);
xlabel('Time (min)');

ff(800,800);
quiver3(d_x(1:end-1),d_y(1:end-1),d_z(1:end-1),d_x(2:end),d_y(2:end),d_z(2:end),0);