files = dir2('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data','-r','*.csv');
minFileSize = 1024 * 50; % MB
files = files([files.bytes] >= minFileSize);

