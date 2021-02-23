function [T,W_z] = detect_sleepWake2(T,n)

W = zeros(size(T.odba_max));
for ii = 1:n
    W = W + smoothdata(T.odba_max,'loess',1440/ii);
end
W_z = (W-mean(W))./std(W);

W_bin = zeros(numel(T.odba),1);
W_bin(W_z > 0) = 1;
T.awake = logical(W_bin);
T.asleep = ~W_bin;