% From Webster:
% The optimal algorithm reached after analysis of the 17 records was:
% D = .025 x [(.15T(i - 4) + .15T(i - 3) + .15T(i - 2) + .08T(i - 1) + .21T(i) + .12T(i + 1) + .13T(i + 2)]
% where T(i) represents the maximal epoch value in minute i, etc. If D >= 1.0, the minute was scored wake;
% otherwise it was scored sleep.
function [T,W] = detect_sleepWake(T,varargin)
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
refVals = T.odba_max(useIds);
T.odba_z = (T.odba_max - mean(refVals)) ./ std(refVals);

if nargin == 1 % no varargin
    nFilt = 2;
    nMinutes = 1;
else
    nFilt = varargin{1};
    nMinutes = varargin{2};
end
% maxFilt = medfilt1(T.odba_z,nFilt);
% maxFilt = smoothdata(T.odba_z,'gaussian',nFilt);
maxFilt = T.odba_z;
% trick to change minutes, but removes data instead of combining it?
maxFilt = equalVectors(maxFilt,numel(maxFilt)/nMinutes);

nR = 4:numel(maxFilt)-1;
im = zeros(numel(nR),5);

% default values on top
if nargin < 3
    im(:,1) = 0.15 * maxFilt(nR-3);
    im(:,2) = 0.15 * maxFilt(nR-2);
    im(:,3) = 0.05 * maxFilt(nR-1);
    im(:,4) = 0.20 * maxFilt(nR);
    im(:,5) = 0.20 * maxFilt(nR+1);
    useThresh = 0;
else % user supplied
    im(:,1) = varargin{3} * maxFilt(nR-3);
    im(:,2) = varargin{4} * maxFilt(nR-2);
    im(:,3) = varargin{5} * maxFilt(nR-1);
    im(:,4) = varargin{6} * maxFilt(nR);
    im(:,5) = varargin{7} * maxFilt(nR+1);
    useThresh = varargin{8};
end

W = [NaN(3,1);sum(im,2);NaN(1,1)]; % multiplier determines threshold, determined emperically
W = interp1(1:numel(maxFilt),W,linspace(1,numel(maxFilt),size(T,1)),'linear')'; % trick to change minutes
W = smoothdata(W,'gaussian',nFilt);
% if numel(W) == 1442 % odd numbers
%     W = W(2:end-1);
% end
Db = zeros(numel(T.odba),1);
Db(W > useThresh) = 1;
% Db = circshift(Db,4);
% Db(1:4) = Db(5);
% Db(end:end-1) = Db(end-2);
T.awake = logical(Db);
T.asleep = ~Db;