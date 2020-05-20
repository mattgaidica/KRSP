function t = smoothTemp(data,paintTemp)
t = data;
for ii = 1:numel(paintTemp)
    t(data == paintTemp(ii)) = NaN;
end

t = inpaint_nans(t);