function colors = mycmap(filename,varargin)
% read image
A = imread(filename);
% set number of elements to return
n = size(A,2);
if ~isempty(varargin)
    n = varargin{1};
end
% create colors (cmap)
useIdx = round(linspace(1,size(A,2),n)); % include end colors
colors = double(squeeze(A(1,useIdx,:))) ./ 255;
% colors = double(squeeze(imresize(A,[1,n]))) ./ 255; % interp end colors