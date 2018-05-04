function prob = get_prob(dists, gammas)
ncol = size(dists, 2);
prob = exp(-gammas * dists);
prob_max = sum(prob, 2);
prob = prob ./ repmat(prob_max, [1, ncol]);
end
