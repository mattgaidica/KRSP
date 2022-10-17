sqkey = readtable('sqkey.txt');
loadPath = '/Users/matt/Dropbox (University of Michigan)/from_box/KRSP Axy Data/Temp';
exportPath = '/Users/matt/Dropbox (University of Michigan)/Biologging/Data';
Tss = makeTss(2014:2020);
%%
AxyDataKey = table;
tableCount = 0;
for iSq = 1:size(sqkey,1)
    T = loadTStruct(iSq,sqkey,Tss);
    if isempty(T)
        continue;
    end
    tableCount = tableCount + 1;
    newFilename = sqkey.filename{iSq};
    newFilename = [newFilename(1:end-3) 'csv'];
    
    try
        T.odba_max = [];
    catch
    end
    try
        T.odba_mean = [];
    catch
    end
    try
        T.odba_z = [];
    catch
    end
    writetable(T,fullfile(exportPath,newFilename));
    
    AxyDataKey.squirrel_id{tableCount} = sqkey.squirrel_id(iSq);
    AxyDataKey.filename{tableCount} = newFilename;
    AxyDataKey.axyCollar{tableCount} = sqkey.axy{iSq};
    AxyDataKey.tag_left{tableCount} = sqkey.tag_left{iSq};
    AxyDataKey.tag_right{tableCount} = sqkey.tag_right{iSq};
    AxyDataKey.grid{tableCount} = sqkey.grid{iSq};
    AxyDataKey.midden{tableCount} = sqkey.midden{iSq};
    AxyDataKey.season{tableCount} = sqkey.season{iSq};
    AxyDataKey.sex{tableCount} = sqkey.sex{iSq};
    AxyDataKey.sex_status{tableCount} = sqkey.sex_status{iSq};
    AxyDataKey.treatment{tableCount} = sqkey.treatment{iSq};
    AxyDataKey.year{tableCount} = sqkey.year(iSq);
    AxyDataKey.experimenter{tableCount} = sqkey.source{iSq};
    AxyDataKey.rec_minutes{tableCount} = sqkey.rec_mins(iSq);
end

% write benKey
writetable(AxyDataKey,fullfile(exportPath,'AxyDataKey_20220730.csv'));