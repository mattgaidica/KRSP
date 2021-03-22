if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';

    odba_nest = [];
    odba_out = [];
    odba_nest_norm = [];
    odba_out_norm = [];
    xs = [];
    female = [];
    sqCount = 0;
    doyCount = 0;
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            if isValidT(T,true)
                sqCount = sqCount + 1;
                dtdoys = day(T.datetime,'dayofyear');
                for doy = unique(dtdoys)'
                    doyCount = doyCount + 1;
                    nest_samples = T.odba(strcmp(T.nest,'Nest') & dtdoys==doy);
                    odba_nest(doyCount) = mean(nest_samples);
                    odba_nest_norm(doyCount) = sum(nest_samples) / numel(nest_samples);
                    out_samples = T.odba(strcmp(T.nest,'Out') & dtdoys==doy);
                    odba_out(doyCount) = mean(out_samples);
                    odba_out_norm(doyCount) = sum(out_samples) / numel(out_samples);
                    xs(doyCount) = doy;
                    female(doyCount) = strcmp(sqkey.sex{iSq},'F');
                end
            end
        end
    end
    do = false;
end

colors = lines(3);
ms = 10;
close all
ff(1200,800,2);
subplot(211);
plot(xs,odba_out,'o','color','k','markersize',ms);
hold on;
plot(xs,odba_nest,'x','color','r','markersize',ms);
plot(xs(female==1),odba_out(female==1),'m.');
xlim([1 366]);
xlabel('doy');
ylabel('odba');
title('mean');

subplot(212);
plot(xs,odba_out_norm,'o','color','k','markersize',ms);
hold on;
plot(xs,odba_nest_norm,'x','color','r','markersize',ms);
plot(xs(female==1),odba_out_norm(female==1),'m.');
xlim([1 366]);
xlabel('doy');
title('normalized');
ylabel('odba');