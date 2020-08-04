if doOnce % at open
    sq_em = readtable('sqkey_emily.txt');
    sq_sa = readtable('/Users/matt/Documents/Data/KRSP/sqkey_sarah_AllAxyMetaData.csv');
    sq_bd = readtable('/Users/matt/Documents/MATLAB/KRSP/2019modexport.csv');
    sq_kr = readtable('/Volumes/Seagate Expansion Drive/2020 XZY data/2020 AXY logsheet.csv');
    
    h_em = sq_em.Properties.VariableNames;
    for ii = 1:numel(h_em)
        sq_em.Properties.VariableNames{ii} = lower(h_em{ii});
    end
    
    h_sa = sq_sa.Properties.VariableNames;
    for ii = 1:numel(h_sa)
        sq_sa.Properties.VariableNames{ii} = lower(h_sa{ii});
    end
    
    h_bd = sq_bd.Properties.VariableNames;
    for ii = 1:numel(h_bd)
        sq_bd.Properties.VariableNames{ii} = lower(h_bd{ii});
    end
    
    h_kr = sq_kr.Properties.VariableNames;
    for ii = 1:numel(h_kr)
        sq_kr.Properties.VariableNames{ii} = lower(h_kr{ii});
    end
    
    % for debugging, used on studd and westrick
% % % %     h_em = sq_em.Properties.VariableNames;
% % % %     h_sa = sq_sa.Properties.VariableNames;
% % % %     h = unique([h_em(:)' h_sa(:)']);
% % % %     for ii = 1:numel(h)
% % % %         in_em = ismember(h{ii},h_em);
% % % %         in_sa = ismember(h{ii},h_sa);
% % % %         fprintf('em%i sa%i %i %s\n',in_em,in_sa,in_em&in_sa,h{ii});
% % % %     end
    
    doOnce = false;
end

moveFrom = '/Users/matt/Box Sync/KRSP Axy Data/Temp/_misfit';
moveTo = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
movefiles = dir(fullfile(moveFrom,'*.mat'));
for iFile = 1:numel(movefiles)
    movefile(fullfile(moveFrom,movefiles(iFile).name),fullfile(moveTo,movefiles(iFile).name));
end

x = {}; % axy for BD data
for ii = 1:numel(sq_bd.axy_)
    x{ii} = num2str(sq_bd.axy_(ii));
end

% fill KR 
xsex = {};
xid = [];
for ii = 1:numel(sq_kr.axy_)
    xsex{ii,1} = 'F'; % all females
    xid(ii,1) = NaN; % ids dont exist yet
end

sqkey = table;
sqkey.axy = [sq_em.axy;sq_sa.axy;x';sq_kr.axy_];
sqkey.filename = [sq_em.filename;sq_sa.rawfilename;sq_bd.filename;sq_kr.datafilename];
sqkey.grid = [sq_em.grid;sq_sa.grid;sq_bd.grid;sq_kr.grid];
sqkey.midden = [sq_em.midden;sq_sa.midden_collar;sq_bd.midden;sq_kr.midden];
sqkey.season = [sq_em.season;sq_sa.deploymentperiod;sq_bd.axyperiod_beginningorlac_;sq_kr.axyperiod_beginningorlac_];
sqkey.sex = [sq_em.sex;sq_sa.sex;sq_bd.partdate_litter1_;xsex];
sqkey.sex_status = strings(size(sqkey,1),1);
sqkey.squirrel_id = [sq_em.squirrel_id;sq_sa.sqid;sq_bd.squirrelid;xid];
sqkey.tag_left = [sq_em.taglft;sq_sa.taglt;sq_bd.momtagl;sq_kr.momtagl];
sqkey.tag_right = [sq_em.tagrt;sq_sa.tagrt;sq_bd.momtagr;sq_kr.momtagr];
sqkey.treatment = [sq_em.treatment;sq_sa.treatment;sq_bd.playbacktreatment;sq_kr.playbacktreatment];
sqkey.year = [sq_em.year;sq_sa.year;year(sq_bd.datestart);repmat(2020,[size(xsex,1),1])];
sqkey.source = [repmat({'ES'},size(sq_em.year));repmat({'SW'},size(sq_sa.year));...
    repmat({'BD'},size(sq_bd.datestart));repmat({'KR'},size(xsex))];

% fix sa filenames
lookPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
lookFiles = dir(fullfile(lookPath,'*.mat'));
count = 0;

for jj = 1:size(sq_sa,1)
    if numel(sq_sa.rawfilename{jj}) >= 5
        fileId = find(contains({lookFiles.name},sq_sa.rawfilename{jj}(1:end-4)));
        if ~isempty(fileId) && contains(sqkey.filename{jj},num2str(sqkey.year(jj)))
            sqkey.filename{strcmp(sqkey.filename,sq_sa.rawfilename{jj})} = lookFiles(fileId).name;
            count = count + 1;
        else
            % try strrep
            fn = strrep(sq_sa.rawfilename{jj}(1:end-4),'ctrl','control');
            fn = strrep(fn,'lactation','lac');
            fileId = find(contains({lookFiles.name},fn));
            if ~isempty(fileId) && contains(sqkey.filename{jj},num2str(sqkey.year(jj)))
                sqkey.filename{strcmp(sqkey.filename,sq_sa.rawfilename{jj})} = lookFiles(fileId).name;
                count = count + 1;
            else
                if sq_sa.year(jj) == 2017
                    altFilestr = [sq_sa.grid{jj},'-',sq_sa.midden_collar{jj},'-',num2str(sq_sa.sqid(jj))];
                    fileId = find(contains({lookFiles.name},altFilestr) & contains({lookFiles.name},'2017'));
                    if ~isempty(fileId)
                        sqkey.filename{strcmp(sqkey.filename,sq_sa.rawfilename{jj})} = lookFiles(fileId).name;
                        count = count + 1;
                    end
                end
            end
        end
    end
end

for jj = 1:size(sq_bd,1)
    if numel(sq_bd.filename{jj}) >= 5
        fileId = find(contains({lookFiles.name},sq_bd.filename{jj}(1:end-4)));
        if ~isempty(fileId)
            sqkey.filename{strcmp(sqkey.filename,sq_bd.filename{jj})} = lookFiles(fileId).name;
        end
    end
end

for jj = 1:size(sq_kr,1)
    if numel(sq_kr.datafilename{jj}) >= 5
        fileId = find(contains({lookFiles.name},sq_kr.datafilename{jj}(1:end-4)));
        if ~isempty(fileId)
            sqkey.filename{strcmp(sqkey.filename,sq_kr.datafilename{jj})} = lookFiles(fileId).name;
        end
    end
end

% All: columns to lower()?

% manual override
sqkey.filename{strcmp(sqkey.filename,'JO-A4-13228-ctrl-lactation.csv')} =... 
    'JO-A4-13228-control-preg-lac-classified_nest__20160404.mat';
% % % % sqkey.filename{strcmp(sqkey.filename,'JO-A4-13228-ctrl-preg.csv')} =...
% % % %     'JO-A4-13228-control-preg-preg-classified_nest__20160314.mat';
sqkey.filename{strcmp(sqkey.filename,'JO-01-21414-GC-preg 2.csv')} =...
    'JO-01-21414-GC-preg-classified_nest__20160404.mat';
sqkey.filename{strcmp(sqkey.filename,'JO-B.5.-13229-ctrl-lactation.csv')} =...
    'JO-B.5.-13229-pregcontrol-lac-classified_nest__20160501.mat';
sqkey.filename{strcmp(sqkey.filename,'JO-B.5.-13229-ctrl-preg.csv')} =...
    'JO-B.5.-13229-pregcontrol-preg-classified_nest__20160410.mat';
sqkey.filename{strcmp(sqkey.filename,'JO-B.12-12700-GC-lactation.csv')} =...
    'JO-B.12-12700-none-lac-classified_nest__20160427.mat';
% % % % posIds = find(strcmp(sqkey.filename,'JO-F.18.-20281-control-lac-classified__20170517.mat'));
% % % % sqkey.filename{posIds(1)} = 'JO-F.18.-20281-pregcontrol-lac-classified_nest__20160417.mat';
% sqkey.filename{strcmp(sqkey.filename,'BJD15_1.csv')} =...
%     'JO-F.18.-20281-control-lac-classified__20170517.mat';

% now clean up
for ii = 1:numel(sqkey.filename)
    if numel(sqkey.filename{ii}) < 5 || ~strcmp(sqkey.filename{ii}(end-3:end),'.mat')
        sqkey.filename{ii} = '';
    end
end
% move misfit files without row in key
mkPath = fullfile(lookPath,'_misfit');
if ~isfolder(mkPath)
    mkdir(mkPath);
end
for ii = 1:numel(lookFiles)
    if ~any(contains(sqkey.filename,lookFiles(ii).name))
        movefile(fullfile(lookPath,lookFiles(ii).name),fullfile(mkPath,lookFiles(ii).name))
        disp(['moving ',lookFiles(ii).name]);
    end
end

% make sex_status, convert sex to M or F
% no entires for breeding right now, leaving for future
statusStr = {'','breeding','nbreeding','pregnant','lactating'}; % verb
for iSq = 1:size(sqkey,1)
    sqkey.sex_status{iSq} = statusStr{1}; % init
    if strcmp(sqkey.sex{iSq},'F-Lac') % Emily's data only
        sqkey.sex{iSq} = 'F';
        sqkey.sex_status{iSq} = statusStr{5};
    end
    if any(strcmp(sqkey.season{iSq},{'lactation','LL','LAC'})) % Emily's data only
        sqkey.sex_status{iSq} = statusStr{5};
    end
    if any(strcmp(sqkey.season{iSq},{'PREG/LAC','pregnancy'})) % Emily's data only
        sqkey.sex_status{iSq} = statusStr{4};
    end
    if any(strcmp(sqkey.season{iSq},{'non-breeding'})) % Emily's data only
        sqkey.sex_status{iSq} = statusStr{3};
    end
end

% just do BD2019 data separately
for iSq = 1:size(sqkey,1)
    if strcmp(sqkey.source{iSq},'BD')
        if strcmp(sqkey.sex{iSq},'Male')
            sqkey.sex{iSq} = 'M';
        else
            sqkey.sex_status{iSq} = sqkey.sex{iSq};
            sqkey.sex{iSq} = 'F';
        end
    end
end

writetable(sqkey,'sqkey');