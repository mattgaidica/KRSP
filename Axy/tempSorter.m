if do
    T = readtable('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/Emily Studds Axy/April 2015 - Week 1/B10_Apr23_2015.csv');
%     T = readtable('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/Emily Studds Axy/Early December 2015/E5_Dec1_2015.csv');
%     T = readtable('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/Emily Studds Axy/Late November 2016/G6_Nov25_2016.csv');
    do = 0;
end
perieventNest = removeShortTransitions(nest,60*15); % optional
tempSorter_init(temp,perieventNest,30,false);
% saveas(gcf,fullfile(savePath,'tempSorter-selection-tool.jpg'));