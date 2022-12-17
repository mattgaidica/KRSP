function temp = filterTemp(T,nSmooth)
% find repeating values and fill them in
temp = T.temp;
[C,ia,ic] = unique(temp);
freq = accumarray(ic,1) / numel(temp);
repeatedValues = find(freq > 0.05); % 10%
for ii=1:numel(repeatedValues)
    temp(temp==C(repeatedValues(ii))) = NaN;
end
temp = inpaint_nans(temp,4); % spring method
temp = smoothdata(temp,'gaussian',nSmooth);
