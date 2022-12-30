function [binNestSense,sense] = nestSenseAlg(temp,odba,nest,senseParams)
sense = {};

% senseParams.thresh = ;
% senseParams.w_kmeans = ;
% senseParams.w_temp = ;
% senseParams.w_tempGrad = ;
% senseParams.w_odba = ;
% senseParams.sm_tempGrad = ;
% senseParams.sm_odba = ;
% senseParams.offset_odba = ;

temp = temp-smoothdata(temp,'movmean',86400 * 3); % subtract moving mean for large changes over time
sense.tempZ = senseParams.w_temp*normalize(temp);

sense.tempGrad = gradient(smoothdata(temp,'gaussian',senseParams.sm_tempGrad * 60));
sense.tempGrad = senseParams.w_tempGrad * normalize(sense.tempGrad .* abs(sense.tempGrad)); % square it for some drama:

sense.smOdba = smoothdata(normalize(odba) + senseParams.offset_odba,'gaussian',senseParams.sm_odba * 60); % normal inside smooth
sense.invOdba = senseParams.w_odba * -sense.smOdba;

% only using this means temp grad and ODBA *together* are important
sense.tempGradObda = sense.tempGrad .* sense.smOdba; % this does NOT get normalized, !!unweighted?

% ODBA should be more important around arbitrary temperatures
sense.tempRange = normalize(-abs(mean(temp)-temp),'range',[0,1]) .* sense.invOdba;

% bias towards original classification
sense.kmeansNest = senseParams.w_kmeans * normalize(nest,'range',[-1 1]);

sense.nest = sense.tempZ + sense.tempGradObda + sense.invOdba + sense.tempRange + sense.kmeansNest;
sense.nest = smoothdata(sense.nest,'gaussian',60 * 5);
% sense.nest = sense.nest./sqrt(abs(sense.nest)); % redlwuce large z-scores
sense.nest = normalize(sense.nest);

binNestSense = normalize(sense.nest > 0,'range'); % turn into binary in/out nest
binNestSense(abs(sense.nest) < senseParams.thresh) = NaN; % set small values to ambiguous
binNestSense = fillmissing(binNestSense,'previous'); % fill with previous
binNestSense = fillmissing(binNestSense,'next'); % in case beginning is NaN