loadPath = '/Users/matt/Dropbox (University of Michigan)/from_box/KRSP Axy Data/Temp';
sqkey = readtable('sqkey.txt');
debExportPath = '/Users/matt/Documents/MATLAB/KRSP/Deb/export';
if do
    weather = readtable('HainesJunction_DailyTemps_Master.csv');
    Tss = makeTss(2014:2020);
    use_iSq = find(strcmp(sqkey.sex_status,'lactating'));
    for iSq = use_iSq'
        T = loadTStruct(iSq,sqkey,Tss);
        if isempty(T)
            continue;
        end
        renameFile = sprintf('%i_%03d_%s.csv',sqkey.squirrel_id(iSq),iSq,datestr(sqkey.rec_startDate(1),'yyyymmDD'));
        writetable(T,fullfile(debExportPath,renameFile));
        disp(iSq);
    end
    do = false;
end