% dataFile = '/Users/matt/Documents/Data/KRSP/doi_10.5061_dryad.1s1m8r7__v1/axy_data_for_ODBA_window_size.csv';
dataFile = '/Users/matt/Documents/Data/KRSP/doi_10.5061_dryad.1s1m8r7__v1/squirrelAxy_2s.csv';
A = readtable(dataFile);

curDate = A.date(1);
disp(curDate);
for ii = 1:size(A,1)
    if curDate ~= A.date(ii)
        curDate = A.date(ii);
        disp(curDate);
    end
end