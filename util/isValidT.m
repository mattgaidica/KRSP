function isValid = isValidT(T,doTemp)
doDebug = true;

isValid = true;

if std(T.odba) < 0.1
    isValid = false;
    if doDebug
        disp('Failed: Std of ODBA too low');
    end
end

% !! Hmm these all depend on Nest?
if doTemp % all Nest depends on temp
    % out ODBA is above a reasonable threshold
    if mean(T.odba(strcmp(T.nest,'Out'))) < 0.15
        isValid = false;
        if doDebug
            disp('Failed: Out of nest ODBA low');
        end
    end
    % deprecated this because it fails for preg females, which may have
    % lots of nest movement?
    % out ODBA is not significantly greater than in ODBA
% % % %     if mean(T.odba(strcmp(T.nest,'Out')))/2 < mean(T.odba(strcmp(T.nest,'Nest')))
% % % %         isValid = false;
% % % %         if doDebug
% % % %             disp('Failed: Nest ODBA suspicious');
% % % %         end
% % % %     end
    % not mostly or hardly Nest categories
    if sum(strcmp(T.nest,'Nest'))/numel(T.nest) > 0.9 || sum(strcmp(T.nest,'Nest'))/numel(T.nest) < 0.1
        isValid = false;
        if doDebug
            disp('Failed: poor Nest distribution');
        end
    end
    % temperature exists
    if std(T.temp) == 0
        isValid = false;
        if doDebug
            disp('Failed: std(temp) = 0');
        end
    end
end