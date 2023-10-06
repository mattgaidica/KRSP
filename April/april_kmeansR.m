function r = april_kmeansR(a)
% Define the number of clusters
b = 2;

% Perform kmeans clustering
[~,C] = kmeans(a, b);

% Calculate the midpoint between the two cluster centers
r = (C(1,:) + C(2,:)) / 2;