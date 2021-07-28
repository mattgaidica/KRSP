if do
    T = readtable('nestAsleepOverlap_v2.csv');
end
clc
lme = fitlme(T,'in_asleep ~is_female+meanSeason+(1|squirrelId)')