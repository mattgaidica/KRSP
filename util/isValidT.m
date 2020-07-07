function isValid = isValidT(T,doTemp)
isValid = true;
% not mostly or hardly Nest categories
if sum(strcmp(T.nest,'Nest'))/numel(T.nest) > 0.8 || sum(strcmp(T.nest,'Nest'))/numel(T.nest) < 0.2
    isValid = false;
end

% out ODBA is not significantly greater than in ODBA
if mean(T.odba(strcmp(T.nest,'Out')))/2 < mean(T.odba(strcmp(T.nest,'Nest')))
    isValid = false;
end

% out ODBA is above a reasonable threshold
if mean(T.odba(strcmp(T.nest,'Out'))) < 0.15
    isValid = false;
end

% out ODBA is above a reasonable threshold
if std(T.temp) == 0 && doTemp
    isValid = false;
end