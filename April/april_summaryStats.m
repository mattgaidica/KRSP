function axy = april_summaryStats(axy)

% Compute derived variables
axy.Xd = axy.x - axy.Xs;
axy.Yd = axy.y - axy.Ys;
axy.Zd = axy.z - axy.Zs;
axy.odba = abs(axy.Xd) + abs(axy.Yd) + abs(axy.Zd);

% Compute Dx, Dy, Dz using the diffX function
n = height(axy);
Dx = nan(n, 1); Dy = Dx; Dz = Dx; % Preallocate
for i = 1:(n-10)
    Dx(i) = diffX(axy.Xd(i:i+10));
    Dy(i) = diffX(axy.Yd(i:i+10));
    Dz(i) = diffX(axy.Zd(i:i+10));
end
axy.Dx = Dx;
axy.Dy = Dy;
axy.Dz = Dz;

% Compute Dodba
axy.Dodba = axy.Dx + axy.Dy + axy.Dz;
axy.Dodba = fillmissing(axy.Dodba, 'linear', 'SamplePoints', 1:height(axy));

% Compute maxY and Todba
axy.maxY = movmax(abs(axy.Yd), [4 0], 'Endpoints', 'fill', 'SamplePoints', 1:height(axy));
axy.Todba = movsum(axy.odba, [10 0], 'Endpoints', 'fill', 'SamplePoints', 1:height(axy));
axy.maxY = fillmissing(axy.maxY, 'linear', 'SamplePoints', 1:height(axy));
axy.Todba = fillmissing(axy.Todba, 'linear', 'SamplePoints', 1:height(axy));

end

function r = diffX(data)
    ab = zeros(1, 11);
    for i = 2:11
        ab(i) = data(i) - data(i-1);
    end
    r = sum(abs(ab));
end

