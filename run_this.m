close all
ff(1200,400,2);
plot(T.odba);
hold on;
T = detect_sleepWake(T,2);
plot(T.awake,'r-')
plot(strcmp(T.nest,'Out')+1);
ylabel('ODBA');
xlabel('time (s)');
legend('ODBA','awake = 1','out of nest = 1');
xlim([1 size(T,1)]);
set(gca,'fontsize',14);