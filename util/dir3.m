function T = dir3(dn,varargin)
inputVals = {'',true}; % filetype, recursive
inputVals(1:nargin-1) = varargin;

if inputVals{2}
    dout = dir2(dn,inputVals{1},'-r');
else
    dout = dir2(dn,inputVals{1});
end

T = struct2table(dout);
T.name = string(T.name);
T.folder = string(T.folder);
T.date = datetime(T.date);
T.fullfile = T.name; % temporary
rmIds = [];
for ii = 1:size(dout,1)
    [~,name,ext] = fileparts(dout(ii).name);
    filename = [name,ext];
    T.name(ii) = filename;
    T.fullfile(ii) = fullfile(T.folder(ii),T.name(ii));
    if any(strcmp(filename(1),{'.','$'})) % rm hidden/recycle
        rmIds = [rmIds ii]; %#ok<AGROW> 
    end
end
T(rmIds,:) = [];