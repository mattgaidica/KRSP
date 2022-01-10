%% how many records have files? This is really the starting point
iFile = 0;
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq})
        iFile = iFile + 1;
    end
end
fprintf("%i no file\n",iFile);

%%
longevity = readtable('krsp_longevity.csv');
cone_counts = readtable('krsp_cone_counts.csv');
midden_cones = readtable('krsp_midden_cones.csv');
iValid = 0;
iLon = 0;
iMidden = 0;
iCones = 0;
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq})
        continue;
    end
%     T = loadTStruct(iSq,sqkey,Tss);
    if sqkey.isValid(iSq)
        iValid = iValid + 1;
        lonId = find(longevity.squirrel_id == sqkey.squirrel_id(iSq));
        if ~isempty(lonId)
            iLon = iLon + 1;
        end
        midId = find(midden_cones.squirrel_id == sqkey.squirrel_id(iSq) &...
            midden_cones.year == sqkey.year(iSq));
        if ~isempty(midId)
            iMidden = iMidden + 1;
        end
    end
end
fprintf("%i valid\n",iValid);
fprintf("%i longevity\n",iLon);
fprintf("%i middens\n",iMidden);

%%
figure;
useId = find(RITable.season == 4 & RITable.is_mast == 0);
x = RITable.age(useId);
y = RITable.RI(useId);
[fitresult, gof] = poly2fit(x,y);
yyaxis left;
plot(fitresult,'k-',x,y,'ko');
hold on;
scatter(x,y,'k','filled');
set(gca,'ycolor','k');

yyaxis right;
useId = find(RITable.season == 4 & RITable.is_mast == 1);
x = RITable.age(useId);
y = RITable.RI(useId);
[fitresult, gof] = poly2fit(x,y);
plot( fitresult,'r-',x,y,'ro');
hold on;
scatter(x,y,'r','filled');
set(gca,'ycolor','r');