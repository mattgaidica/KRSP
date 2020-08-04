function cmap = magma(varargin)
cmapPath = '/Users/matt/Documents/MATLAB/KRSP/util/magma_trans.png';
if ~isempty(varargin) % set nColors
    cmap = mycmap(cmapPath,varargin{1});
else
    cmap = mycmap(cmapPath);
end