sqkey = readtable('sqkey');
filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
saveDir = fullfile(filePath,'_figs');
if ~isfolder(saveDir)
    mkdir(saveDir);
end

for iSq = 1:size(sqkey,1)
    if ~isempty(sqkey.filename{iSq})
        load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
        if isValidT(T,false)
            disp(sqkey.filename{iSq});
            h = ff(1200,500);
            plot(T.datetime,T.odba);
            ylabel('ODBA');
            yyaxis right;
            plot(T.datetime,T.temp);
            ylabel('collar temp');
            title({sqkey.filename{iSq},...
                sprintf('ODBA: %1.2f + %1.2f, temp: %1.2f + %1.2f',...
                mean(T.odba),std(T.odba),mean(T.temp),std(T.temp))},'interpreter','none');
            saveas(h,fullfile(saveDir,[sqkey.filename{iSq},'.jpg']));
            close(h);
        end
    end
end