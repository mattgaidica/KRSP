function cmap = jupiter(varargin)
cmapPath = '/Users/matt/Documents/MATLAB/Development/ChoiceTask/LFPs/utils/jupiter.jpg';
if ~isempty(varargin) % set nColors
    cmap = mycmap(cmapPath,varargin{1});
else
    cmap = mycmap(cmapPath);
end