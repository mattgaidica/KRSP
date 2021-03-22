function Tss = makeTss(useYears)
warning ('off','all');
iRow = 1;
Tss = table;
for thisYear = useYears
    [sunrise,sunset,noon] = estSunrise(thisYear);
    
    Tss.sunrise(iRow:iRow+size(sunrise,2)-1) = sunrise';
    Tss.sunset(iRow:iRow+size(sunset,2)-1) = sunset';
    Tss.noon(iRow:iRow+size(noon,2)-1) = noon';
    
    iRow = iRow + size(sunset,2);
end

Tss.day_length = seconds(duration(Tss.sunset-Tss.sunrise));
Tss.year = year(Tss.noon);
Tss.month = month(Tss.noon);
Tss.day = day(Tss.noon);
Tss.doy = day(Tss.noon,'dayofyear');
warning ('on','all');