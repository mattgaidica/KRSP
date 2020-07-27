% The optimal algorithm reached after analysis of the 17 records was:
% D = .025 x [(.15T(i - 4) + .15T(i - 3) + .15T(i - 2) + .08T(i - 1) + .21T(i) + .12T(i + 1) + .13T(i + 2)]
% where T(i) represents the maximal epoch value in minute i, etc. If D >= 1.0, the minute was scored wake;
% otherwise it was scored sleep.
function [T,warn] = detect_sleepWake(T)
% issues with sleep detection
% 1. a common threshold only works on homogenous data
% 2. a binary classifier assumes certainty
% 3. a classifier without history assumes sleep state is indepedent of wake

% get sunrise
sunrise = secDay(sunriseSunset(T.datetime(1))); % just use first day
allSecs = secDay(T.datetime);
% select ~3hrs before sunrise
testSpan = 3*60*60; % 3 hours
sunriseOffset = 15*60; % 15 minutes
useIds = allSecs < (sunrise-sunriseOffset) & allSecs > (sunrise-sunriseOffset-testSpan);
refVals = T.odba(useIds);
T.odba_z = (T.odba - mean(refVals)) ./ std(refVals);

nR = 5:numel(T.odba_z)-2;
nFilt = 2;
maxFilt = medfilt1(T.odba_z,nFilt);
im = zeros(numel(nR),7);
im(:,1) = 0.15 * maxFilt(nR-4);
im(:,2) = 0.15 * maxFilt(nR-3);
im(:,3) = 0.15 * maxFilt(nR-2);
im(:,4) = 0.08 * maxFilt(nR-1);
im(:,5) = 0.21 * maxFilt(nR);
im(:,6) = 0.21 * maxFilt(nR+1);
im(:,7) = 0.21 * maxFilt(nR+2);

W = sum(im,2); % multiplier determines threshold, determined emperically
Db = zeros(numel(T.odba_max),1);
Db(W > 0) = 1;
Db = circshift(Db,4);
Db(1:4) = Db(5);
Db(end:end-1) = Db(end-2);
T.awake = logical(Db);
T.asleep = ~Db;