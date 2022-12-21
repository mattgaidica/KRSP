function [binNestSense,nestSense] = nestSenseAlg(temp,odba,wArr)

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
comp2 = gradient(smoothdata(temp,'gaussian',gradSm));
comp2 = w_tempGrad*normalize(comp2.*abs(comp2));
comp3 = w_odba*-smoothdata(normalize(odba),'gaussian',odbaSm);
nestSense = comp1+comp2+comp3;
nestSense = smoothdata(normalize(nestSense),'gaussian',1);
binNestSense = normalize(nestSense>useThresh,'range');