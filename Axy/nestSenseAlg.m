function [binNestSense,nestSense,alg_temp,alg_tempRange,alg_invOdba,alg_tempGradObda,alg_kmeansNest,smOdba]...
    = nestSenseAlg(temp,odba,nest,wArr,zOffset)
useThresh = 0.2;
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

% temp has to sweep across days (like original k-means alg)
temp = temp-smoothdata(temp,'movmean',86400*3);
alg_temp = w_temp*normalize(temp);

alg_tempGrad = gradient(smoothdata(temp,'gaussian',gradSm*60));
alg_tempGrad = w_tempGrad*normalize(alg_tempGrad.*abs(alg_tempGrad)); % square it for some drama:

smOdba = smoothdata(normalize(odba)+zOffset,'gaussian',odbaSm*60); % normal inside smooth
alg_invOdba = w_odba*-smOdba;

% only using this means temp grad and ODBA *together* are important
alg_tempGradObda = alg_tempGrad.*smOdba; % this does NOT get normalized, !!unweighted?

% ODBA should be more important around arbitrary temperatures
alg_tempRange = normalize(-abs(mean(temp)-temp),'range',[0,1]).*alg_invOdba;

% bias towards original classification
alg_kmeansNest = w_kmeans*normalize(nest,'range',[-1 1]);

nestSense = alg_temp + alg_tempGradObda + alg_invOdba + alg_tempRange + alg_kmeansNest;
nestSense = smoothdata(nestSense,'gaussian',60*5);
% nestSense = nestSense./sqrt(abs(nestSense)); % redlwuce large z-scores
nestSense = normalize(nestSense);

% nestGrad = normalize(smoothdata(gradient(nestSense),'gaussian',60*5));

binNestSense = normalize(nestSense>0,'range'); % turn into binary in/out nest
binNestSense(abs(nestSense) < useThresh) = NaN; % set small values to ambiguous
binNestSense = fillmissing(binNestSense,'previous'); % fill with previous
binNestSense = fillmissing(binNestSense,'next'); % in case beginning is NaN