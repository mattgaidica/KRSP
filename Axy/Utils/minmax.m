function [mmag,mmin,mmax] = minmax(arr)
mmin = min(arr);
mmax = max(arr);
mmag = mmax - mmin;
if round(mmin) == mmin && round(mmax) == mmax
    fprintf("%i-%i, diff=%i\n",mmin,mmax,mmag);
else
    fprintf("%1.3f-%1.3f, diff=%1.3f\n",mmin,mmax,mmag);
end