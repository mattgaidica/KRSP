function [mmag,mmin,mmax] = minmax(arr)
mmin = min(arr);
mmax = max(arr);
mmag = mmax - mmin;