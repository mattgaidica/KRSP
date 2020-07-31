if doOnce % at open
    sq_em = readtable('sqkey_emily.txt');
    sq_sa = readtable('/Users/matt/Documents/Data/KRSP/sqkey_sarah_AllAxyMetaData.csv');
    sq_bd = readtable('/Users/matt/Documents/MATLAB/KRSP/2019modexport.csv');
    
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
    
    % for debugging, used on studd and westrick
    h_em = sq_em.Properties.VariableNames;
    h_sa = sq_sa.Properties.VariableNames;
    h = unique([h_em(:)' h_sa(:)']);
    for ii = 1:numel(h)
        in_em = ismember(h{ii},h_em);
        in_sa = ismember(h{ii},h_sa);
        fprintf('em%i sa%i %i %s\n',in_em,in_sa,in_em&in_sa,h{ii});
    end
    
    doOnce = false;
end

moveFrom = '/Users/matt/Box Sync/KRSP Axy Data/Temp/_misfit';
moveTo = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
movefiles = dir(fullfile(moveFrom,'*.mat'));
for iFile = 1:numel(movefiles)
    movefile(fullfile(moveFrom,movefiles(iFile).name),fullfile(moveTo,movefiles(iFile).name));
end

sqkey = table;
sqkey.axy = [sq_em.axy;sq_sa.axy;cellfun(@num2str,x,'UniformOutput',false)];
sqkey.filename = [sq_em.filename;sq_sa.rawfilename;sq_bd.filename];
sqkey.grid = [sq_em.grid;sq_sa.grid;sq_bd.grid];
sqkey.midden = [sq_em.midden;sq_sa.midden_collar;sq_bd.midden];
sqkey.season = [sq_em.season;sq_sa.deploymentperiod;sq_bd.axyperiod_beginningorlac_];
sqkey.sex = [sq_em.sex;sq_sa.sex;sq_bd.partdate_litter1_];
sqkey.sex_status = strings(size(sqkey,1),1);
sqkey.squirrel_id = [sq_em.squirrel_id;sq_sa.sqid;sq_bd.squirrelid];
sqkey.tag_left = [sq_em.taglft;sq_sa.taglt;sq_bd.momtagl];
sqkey.tag_right = [sq_em.tagrt;sq_sa.tagrt;sq_bd.momtagr];
sqkey.treatment = [sq_em.treatment;sq_sa.treatment;sq_bd.playbacktreatment];
sqkey.year = [sq_em.year;sq_sa.year;year(sq_bd.datestart)];
sqkey.source = [repmat({'ES'},size(sq_em.year));repmat({'SW'},size(sq_sa.year));...
    repmat({'BD'},size(sq_bd.datestart))];

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
posIds = find(strcmp(sqkey.filename,'JO-F.18.-20281-control-lac-classified__20170517.mat'));
sqkey.filename{posIds(1)} = 'JO-F.18.-20281-pregcontrol-lac-classified_nest__20160417.mat';
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
    if any(strcmp(sqkey.season{iSq},{'lactation','LL'})) % Emily's data only
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

% em1 sa1 1 axy
% em1 sa0 0 census_month
% em1 sa0 0 census_year
% em1 sa0 0 cones_tree
% em1 sa0 0 cones_tree_x1
% em0 sa1 0 dateremoved
% em0 sa1 0 datestart
% em0 sa1 0 dateturnedoff
% em0 sa1 0 daysdeployed
% em1 sa0 0 density
% em0 sa1 0 deploydate
% em0 sa1 0 deployday
% em1 sa0 0 deployment_id
% em0 sa1 0 deploymentperiod
% em0 sa1 0 deploytime
% em0 sa1 0 dose
% em1 sa0 0 filename
% em1 sa0 0 frequency
% em1 sa1 1 grid
% em0 sa1 0 litterid
% em1 sa0 0 locx
% em1 sa0 0 locy
% em1 sa0 0 midden
% em1 sa0 0 midden_1
% em0 sa1 0 midden_collar
% em1 sa0 0 middenyear
% em0 sa1 0 notes
% em0 sa1 0 parturitiondate
% em0 sa1 0 pbend
% em0 sa1 0 pbstart
% em1 sa0 0 post_pb
% em1 sa0 0 pre_pb
% em0 sa1 0 rawfilename
% em0 sa1 0 removeday
% em1 sa0 0 row
% em1 sa0 0 season
% em1 sa1 1 sex
% em0 sa1 0 shake1start
% em0 sa1 0 shake1stop
% em0 sa1 0 shake2start
% em0 sa1 0 shake2stop
% em0 sa1 0 sqid
% em1 sa0 0 squirrel
% em1 sa0 0 squirrel_id
% em1 sa0 0 taglft
% em0 sa1 0 taglt
% em1 sa1 1 tagrt
% em0 sa1 0 timeremoved
% em0 sa1 0 timestart
% em0 sa1 0 timeturnedoff
% em1 sa1 1 treatment
% em0 sa1 0 trtstartday
% em0 sa1 0 trtstopday
% em1 sa0 0 var1
% em1 sa1 1 year
