function [binNestSense,nestSense,comp1,comp2,comp3] = nestSenseAlg(temp,odba,wArr)

if isempty(wArr)
    useThresh = 0;
    w_temp = 3;
    w_tempGrad = 1;
    w_odba = 1;
    gradSm = 60*10;
    odbaSm = 60*10;
else
    useThresh = wArr(1);
    w_temp = wArr(2);
    w_tempGrad = wArr(3);
    w_odba = wArr(4);
    gradSm = wArr(5);
    odbaSm = wArr(6);
end
comp1 = w_temp*normalize(temp);
comp2 = gradient(smoothdata(temp,'gaussian',gradSm*60));
comp2 = w_tempGrad*normalize(comp2.*abs(comp2));
comp3 = w_odba*-normalize(smoothdata(odba,'gaussian',odbaSm*60));
% comp4 = w_temp*smoothdata(normalize(temp),'gaussian',60*60*10);
nestSense = comp1+abs(comp2).*comp3;
nestSense = normalize(smoothdata(nestSense,'gaussian',60*5));

% nestGrad = normalize(smoothdata(gradient(nestSense),'gaussian',60*5));

binNestSense = normalize(nestSense>useThresh,'range');
% binNestSense(abs(nestGrad) < 1) = NaN;
% binNestSense = fillmissing(binNestSense,'nearest');