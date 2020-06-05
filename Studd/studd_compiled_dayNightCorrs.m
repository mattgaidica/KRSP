pThreshs = linspace(.01,0.99,10);
tw = 3600/2;
centerTimes = linspace(-5*3600,-tw,16);
% centerTime = round(3 * 3600);
nSmooth = 90;


% create thresholds based on all nights for subject
if true
    subjects = unique(meta.subjects);
    compiled_night_odba = cell(numel(subjects),1);
    compiled_day_odba = cell(numel(subjects),1);
    for iDay = 1:numel(meta.dates)
        useId = find(strcmp(meta.subjects{iDay},subjects));
        odba_night = smooth(meta.odba_night{iDay},nSmooth);
        compiled_night_odba{useId} = [compiled_night_odba{useId};odba_night];
        odba_day = smooth(meta.odba_day{iDay},nSmooth);
        compiled_day_odba{useId} = [compiled_day_odba{useId};odba_day];
    end
    v_nThresh = [];
    v_dThresh = [];
    for iSubject = 1:numel(subjects)
        for iiN = 1:numel(pThreshs)
            v_night = sort(compiled_night_odba{iSubject});
            v_nThresh(iSubject,iiN) = v_night(round(numel(v_night)*pThreshs(iiN)));
            v_day = sort(compiled_day_odba{iSubject});
            v_dThresh(iSubject,iiN) = v_day(round(numel(v_day)*pThreshs(iiN)));
        end
    end
end

p_mat = [];
r_mat = [];
for iTime = 1:numel(centerTimes)
    for iiN = 1:numel(pThreshs)
        for jjD = 1:numel(pThreshs)
            fprintf('iTime %i, iiN %i, jjD %i\n',iTime,iiN,jjD);
            for iDay = 1:numel(meta.dates)
                useId = find(strcmp(meta.subjects{iDay},subjects));

                odba_night = smooth(meta.odba_night{iDay},nSmooth);
                odba_night_counts(iDay) = sum(odba_night(end+centerTimes(iTime)-tw:...
                    end+centerTimes(iTime)+tw-1) > v_nThresh(useId,iiN));

                odba_day = smooth(meta.odba_day{iDay},nSmooth);
                odba_day_counts(iDay) = sum(odba_day > v_dThresh(useId,jjD)) / meta.days_length(iDay);
            end
            [r,p] = corr(odba_day_counts',odba_night_counts');
            r_mat(iiN,jjD,iTime) = r;
            p_mat(iiN,jjD,iTime) = p;
        end
    end
end

close all
ff(1200,800);
hs = [];
hs(1) = subplot(121);
hs1 = slice(r_mat,[],[],1:numel(centerTimes));
set(hs1,'FaceAlpha',0.6);
colormap(jet);
caxis([-0.3 0.3]);
xticks(linspace(min(xlim),max(xlim),nTicks));
xtickangle(90);
xticklabels(compose('%1.2f',linspace(pThreshs(1),pThreshs(end),nTicks)));
yticks(xticks);
yticklabels(xticklabels);
zticks(min(zlim):max(zlim));
zticklabels(compose('%1.2f',centerTimes/3600));
xlabel('night ODBA frac. threshold');
ylabel('day ODBA frac. threshold');
zlabel('hours relative to waking');
view(140,45);
colorbar;
% shading interp;

hs(2) = subplot(122);
hs2 = slice(p_mat,[],[],1:numel(centerTimes));
set(hs2,'FaceAlpha',0.6);
caxis([0 0.05]);
xticks(linspace(min(xlim),max(xlim),nTicks));
xtickangle(90);
xticklabels(compose('%1.2f',linspace(pThreshs(1),pThreshs(end),nTicks)));
yticks(xticks);
yticklabels(xticklabels);
zticks(min(zlim):max(zlim));
zticklabels(compose('%1.2f',centerTimes/3600));
xlabel('night ODBA frac. threshold');
ylabel('day ODBA frac. threshold');
zlabel('hours relative to waking');
view(140,45);
colorbar;
% shading interp;

ff(900,900);
rows = ceil(sqrt(numel(pThreshs)));
cols = rows;
for ii = 1:numel(centerTimes)
    subplot(rows,cols,ii);
    imagesc(squeeze(p_mat(:,:,ii))');
    colormap(jet);
    caxis([0 0.05]);
    set(gca,'ydir','normal');
    cb = colorbar;
    cb.Label.String = 'p';
    xticks(linspace(min(xlim),max(xlim),nTicks));
    xtickangle(90);
    xticklabels(compose('%1.2f',linspace(pThreshs(1),pThreshs(end),nTicks)));
    yticks(xticks);
    yticklabels(xticklabels);
    title(sprintf('t = %1.2f hrs',centerTimes(ii)/3600));
    set(gca,'fontsize',10);
end

ff(900,900);
rows = ceil(sqrt(numel(pThreshs)));
cols = rows;
for ii = 1:numel(centerTimes)
    subplot(rows,cols,ii);
    imagesc(squeeze(r_mat(:,:,ii))');
    colormap(jet);
    caxis([-.3 .3]);
    set(gca,'ydir','normal');
    cb = colorbar;
    cb.Label.String = 'p';
    xticks(linspace(min(xlim),max(xlim),nTicks));
    xtickangle(90);
    xticklabels(compose('%1.2f',linspace(pThreshs(1),pThreshs(end),nTicks)));
    yticks(xticks);
    yticklabels(xticklabels);
    title(sprintf('t = %1.2f hrs',centerTimes(ii)/3600));
    set(gca,'fontsize',10);
end