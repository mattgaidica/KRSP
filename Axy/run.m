% T=readtable('/Volumes/GAIDICASSD/KRSP/Axy Database/AxyDatabase.xlsx') 
jj = 0;
onDates = string;
onTimes = string;
for ii = 271:653
    jj = jj + 1;
    onDates(jj,:) = "";
    onTimes(jj,:) = "";
    [files,n] = dir3('/Volumes/GAIDICASSD/KRSP/Axy Database/Archive/Emily Studds Axy',T.filename{ii});
    if n == 1
        matchIdx = find(strcmp(T_AxyFiles.filename,files.name(1)));
        if numel(matchIdx) == 1
            onDates(jj,:) = string(T_AxyFiles.startDate(matchIdx),'dd/MM/yyyy');
            onTimes(jj,:) = string(T_AxyFiles.startDate(matchIdx),'HH:mm:ss.S');
        end
    end
end