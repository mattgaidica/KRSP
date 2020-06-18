function secs = secDay(t_datetime)
[~,~,~,h,mi,s] = datevec(t_datetime);
secs = h*3600 + mi*60 + s;