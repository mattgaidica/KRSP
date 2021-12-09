function varargout = plotLeftOverRitht(LR,varargin)
% wrapper to plot the left axis of dual-axis plot over the right axis

%% process input
% input string
if strcmpi(LR,'left')
    LR = 'left';
    val = 1;
elseif strcmpi(LR,'right')
    LR = 'right';
    val = 0;
else
    error('plotLeftOverRitht:Input:LR',"The first input must be either 'left' or 'right'.")
end

% input: axis handle?
if isa(varargin{1},'matlab.graphics.axis.Axes')
    ax = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca;
end

%% main function
% activate left axis
yyaxis( ax, LR);
% call plot
lines = plot(   ax, varargin{:});
% set height
setZData(lines,val);
% set sorting order
set(ax, 'SortMethod', 'depth')

%% output
if nargout > 0
    varargout = ax;
end
end