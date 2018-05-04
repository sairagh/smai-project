function gamma = getGamma(dist, excess)
low = 0; high = 100;
[ndata nclust] = size(dist);
excess = excess * nclust;                                          
if excess == nclust
    gamma = low;
    return
end
while high - low > 1e-4
    mid = (low + high) / 2;
    prob = exp(-mid * dist);
    texcess = sum(prob(:)) / ndata;
    if texcess > excess
        low = mid;
    else
        high = mid;
    end
end
gamma = mid;
dif_gamma = abs(texcess - excess);
end
