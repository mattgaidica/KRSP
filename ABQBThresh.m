threshCount = 0;
testedCount = 0;
% allThresh = [];
for iSq = 1:size(sqkey,1)
    [T,sleepThresh] = loadTStruct(iSq,sqkey,Tss);
    if isempty(T)
        continue;
    end
    
    if strcmp(sqkey.sex_status{iSq},'lactating') | strcmp(sqkey.sex_status{iSq},'pregnant') |...
            strcmp(sqkey.sex_status{iSq},'Pre-pregnancy')
        disp('skipping female');
        continue;
    end
    
    testedCount = testedCount + 1;
    
    if ~isnan(sleepThresh)
        threshCount = threshCount + 1;
        allThresh(threshCount) = sleepThresh;
        saveas(gcf,sprintf("/Users/matt/Documents/MATLAB/KRSP/export/_QBDetect/%3d.jpg",iSq));
        close all;
    else
        fprintf("no thresh %i\n",iSq);
    end
end
%%
close all
ff(400,400);
bins = linspace(0,0.5,24);
histogram(allThresh,bins,'facecolor','k');
threshMed = median(allThresh);
hold on;
[h,p] = lillietest(allThresh); %is it normally distrubted?
xline(threshMed,'r-','linewidth',2);
ylim([0 max(ylim)+6]);
text(threshMed,max(ylim)-5,"\leftarrow" + sprintf("median thresh = %1.3f",threshMed),'fontsize',14,'color','r');
set(gca,'fontsize',14);
xlabel('ODBA');
ylabel('Session Count');
title(sprintf("Valid Thresholds\n%i/%i Rec. Sessions\nLillie test p = %1.2e",...
    numel(allThresh),testedCount,p));

saveas(gcf,fullfile(exportPath,"thresholdDistribution.jpg"));