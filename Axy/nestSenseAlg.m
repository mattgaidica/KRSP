function [binNestSense,nestSense,comp1,comp1_3,comp3,comp2_3,comp4] = nestSenseAlg(temp,odba,nest,wArr)
useThresh = 0;
if isempty(wArr)
%     useThresh = 0;
%     w_temp = 3;
%     w_tempGrad = 1;
%     w_odba = 1;
%     gradSm = 60*10;
%     odbaSm = 60*10;
else
    w_kmeans = wArr(1);
    w_temp = wArr(2);
    w_tempGrad = wArr(3);
    w_odba = wArr(4);
    gradSm = wArr(5);
    odbaSm = wArr(6);
end
% adj_odba = (odba - mean(odba(nest==0))) ./ std(odba);

comp1 = w_temp*normalize(temp);
comp2 = gradient(smoothdata(temp,'gaussian',gradSm*60));
comp2 = w_tempGrad*normalize(comp2.*abs(comp2)); % square it for some drama:

% this fixes the nest guess from temp only then minimizes ODBA in nest
% while maximizing ODBA out of nest
% % rmShortNest = removeShortTransitions(nest,60); % optional
% % shiftNest = fixTempDelay(rmShortNest,odba,temp); % re-align nest
% % % most odba in-nest should be Z<0, most out Z>0
addVals = linspace(0,1,100);
% % resArr = NaN(size(addVals));
% % addIdx = 1;
% % for ii = 1:numel(addVals)
% %     newComp = normalize(odba) + addVals(ii);
% %     resArr(ii,1) = sum(sign(newComp(nest==0))==1)/sum(nest==0);
% %     resArr(ii,2) = sum(sign(newComp(nest==1))==-1)/sum(nest==1);
% %     if resArr(ii,1) >= resArr(ii,2)
% %         addIdx = ii;
% %         break;
% %     end
% % end
% !! figure out how to run code above just once? cache addIdx?
addIdx = 33;
% figure;plot(addVals(1:ii),resArr(1:ii,1));hold on;plot(addVals(1:ii),resArr(1:ii,2));grid on;
smOdba = smoothdata(normalize(odba)+addVals(addIdx),'gaussian',odbaSm*60);
comp3 = w_odba*-smOdba; % normal inside smooth

% only using this means temp grad and ODBA *together* are important
comp2_3 = comp2.*smOdba; % this does NOT get normalized, !!unweighted?

% ODBA should be more important around arbitrary temperatures
comp1_3 = normalize(-abs(mean(temp)-temp),'range',[0,1]).*comp3;

% bias towards original classification
comp4 = w_kmeans*normalize(nest,'range',[-1 1]);

nestSense = comp1+comp2_3+comp3+comp1_3+comp4;
nestSense = smoothdata(nestSense,'gaussian',60*5);
nestSense = nestSense./sqrt(abs(nestSense)); % reduce large z-scores
nestSense = normalize(nestSense);

% nestGrad = normalize(smoothdata(gradient(nestSense),'gaussian',60*5));

binNestSense = normalize(nestSense>useThresh,'range');
% binNestSense(abs(nestGrad) < 1) = NaN;
% binNestSense = fillmissing(binNestSense,'nearest');