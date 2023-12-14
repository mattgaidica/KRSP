function axy = april_behavClass2(axy)
% Define new variables based on conditions
axy.Nest = repmat("Out", height(axy), 1);
axy.Nest(axy.temp11 > axy.center) = "Nest";

axy.Move = repmat("Moving", height(axy), 1);
axy.Move(axy.Dodba <= 1.06) = "NotMoving";

% Adjusted Feed classification
axy.Feed = repmat("Feed", height(axy), 1);
axy.Feed(axy.Move == "NotMoving") = "NotMoving";
axy.Feed(axy.Move ~= "NotMoving" & axy.Todba >= 6.2) = "Forage";

% Adjusted Travel classification
axy.Travel = repmat("Travel", height(axy), 1);
axy.Travel(axy.Feed ~= "Forage") = axy.Feed(axy.Feed ~= "Forage");
axy.Travel(axy.Feed == "Forage" & axy.maxY >= 1.15) = "Forage";

% If using a nestCorrection function
% axy.Nest2 = nestCorrection(axy); % Assuming nestCorrection is defined

% Compute the 'All' variable based on conditions
axy.All = repmat("", height(axy), 1);
axy.All(axy.Nest == "Nest" & axy.Move == "Moving") = "NestMove";
axy.All(axy.Nest == "Nest" & axy.Move == "NotMoving") = "NestNotMove";
axy.All(axy.Nest ~= "Nest") = axy.Travel(axy.Nest ~= "Nest");