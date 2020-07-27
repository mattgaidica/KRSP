tW = 60; % min
sweep = round(linspace(tW+1,720,200));
% nSpace = 60*9;
spaces = round(linspace(60*6,60*12,20));
tSeasons = round(linspace(1,366,5));
seasonDoys = circshift(1:366,60);
unSqs = unique(sq_ids);
% % sq_odba_norm = normalize(sq_odba,2,'scale');
sq_odba_asleep = sq_asleep;
sq_odba_awake = sq_odba_z;
% sq_odba_awake(sq_odba_z < 1) = NaN;
sq_odba_awake(sq_odba > 1 | sq_odba < 0.25) = NaN;
if do
    all_rs = [];
    all_ps = [];
    for iSpace = 1:numel(spaces)
        nSpace = spaces(iSpace);
        for iSeason = 1:4
            disp(iSeason);
            for iSw = 1:numel(sweep)
                useIds = ismember(sq_doys,seasonDoys(tSeasons(iSeason):tSeasons(iSeason+1)));
% %                 useIds = sq_doys >= tSeasons(iSeason) & sq_doys < tSeasons(iSeason+1);
                x = [];
                y = [];
                bIds = sweep(iSw)-tW : sweep(iSw);
                fIds = sweep(iSw)+nSpace-tW : sweep(iSw)+nSpace;
                x = mean(sq_odba_asleep(useIds,bIds),2);
                y = nanmean(sq_odba_awake(useIds,fIds),2);
                [r,p] = corr(x,y,'rows','complete');
                if p < 0.05
                    all_rs(iSeason,iSw,iSpace) = r;
                else
                    all_rs(iSeason,iSw,iSpace) = NaN;
                end
                all_rs(iSeason,iSw,iSpace) = r;
                all_ps(iSeason,iSw,iSpace) = p;
            end
        end
    end
    % do = false;
end

%%
titleLabels = {'winter','spring','summer','fall'};
ff(400,600,2);
for iSeason = 1:4
    subplot(4,1,iSeason);
    theseRs = squeeze(all_rs(iSeason,:,:));
    imagesc(theseRs');
    colorbar;
    caxis([-0.25 0.25]);
    colormap(jupiter);
    yticklabels(round(spaces(yticks)/60));
    xlim([1 size(all_rs,2)]);
    xticks(xlim);
    xticklabels({sprintf('-%ihrs',nSpace/60),...
        sprintf('sunrise',nSpace/60)});
    ylabel('t+hours');
    xlabel('t');
    title(titleLabels{iSeason});
end
% pos corr => high~high, low~low % not expected
% neg corr => high~low, low~high

%% plot
% % % % ff(400,600,2);
% % % % for iSeason = 1:4
% % % %     subplot(4,1,iSeason);
% % % %     theseRs = all_rs(iSeason,:);
% % % %     thesePs = all_ps(iSeason,:);
% % % %     plot(theseRs,'k-','linewidth',2);
% % % %     hold on;
% % % %     plot(find(thesePs < 0.05),theseRs(find(thesePs < 0.05)),'r.');
% % % %     %     plot(nanmean(theseVals) - nanstd(theseVals),'k--');
% % % %     %     plot(nanmean(theseVals) + nanstd(theseVals),'k--');
% % % %     plot(xlim,[0 0],'k:');
% % % %     ylim([-.4 .4]);
% % % %     yticks([min(ylim) 0 max(ylim)]);
% % % %     ylabel('corr coeff');
% % % %     title(sprintf('Q%i %imin win, %ihrs forward corr',iSeason,tW,nSpace/60));
% % % %     xlim([1 size(all_rs,2)]);
% % % %     xticks(xlim);
% % % %     xticklabels({sprintf('-%ihrs',nSpace/60),...
% % % %         sprintf('sunrise',nSpace/60)});
% % % % end