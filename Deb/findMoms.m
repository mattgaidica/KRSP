Litter = readtable(fullfile('R','krsp_litter.csv'));
Juvenile = readtable(fullfile('R','krsp_juvenile.csv'));
Squirrel = readtable(fullfile('R','krsp_squirrel.csv'));
sqkey = readtable('sqkey.txt');
%%
momCount = 0;
juvAge_atStartRec = [];
for iSq = 1:size(sqkey,1)
    litterIds = find(Litter.squirrel_id == sqkey.squirrel_id(iSq));
    for litterId = litterIds'
        litterDate = Litter.fieldBDate(litterId);
        recStart = sqkey.rec_startDate(iSq);
        if ~isnat(litterDate) & ~isnat(recStart)
            momCount = momCount + 1;
            juvAge_atStartRec(momCount) = days(litterDate - recStart);
        end
    end
end
close all
ff(400,400);
histogram(juvAge_atStartRec(juvAge_atStartRec > 0),linspace(0,40,9),'FaceColor','k');
title(sprintf('Litters ≤ 40 days old w/ Axy (n = %i)',sum(juvAge_atStartRec > 0 & juvAge_atStartRec <= 40)));
set(gca,'fontsize',14);
ylabel('# Litters');
xlabel('Age of Litter at Axy Start (days)');
ylim(ylim.*[0 1.1]);
saveas(gcf,'AxyLitters.png');