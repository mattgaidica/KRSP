function [T,n] = dir3(dn,varargin)
inputVals = {'',true,{}}; % filter ("*.csv"), do recursive, rm if contains
inputVals(1:nargin-1) = varargin;
T = table;
n = 0;
if ~ismissing(string(dn)) && ~ismissing(string(inputVals{1}))
    if inputVals{2}
        dout = dir2(char(dn),char(inputVals{1}),'-r');
    else
        dout = dir2(char(dn),char(inputVals{1}));
    end
    
    T = struct2table(dout);
    if ~isempty(dout)
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
            if any(strcmp(filename(1),{'.','$'})) || contains(lower(filename),lower(inputVals{3})) % rm hidden/recycle
                rmIds = [rmIds ii]; %#ok<AGROW> 
            end
        end
        T(rmIds,:) = [];
    end
    n = size(T,1);
end