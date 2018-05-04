function dists = dist_p2c(ktr, kte, clustidx)
nclust = max(clustidx);
nte = size(kte, 1);                         
dists = zeros(nte, nclust);                 
for iclust = 1 : nclust
    ind = (clustidx == iclust);
    card = sum(ind);                       
    dists(:, iclust) = sum(sum(ktr(ind, ind))) / (card^2) - 2 * sum(kte(:, ind), 2) / card;
end
dist_min = min(dists, [], 2);
dists = dists - repmat(dist_min, [1, nclust]);
end
