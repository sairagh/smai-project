function s = apply_soft(model, Kts, prob_ts)
[ntest, nsv, nkernel] = size(Kts);
betas = model.betas;
prob_sv = model.prob_sv;
kaux = zeros(ntest, nsv);
nclust = size(prob_sv, 2);                                                   
for j = 1 : nclust
    kaux = kaux + (prob_ts(:, j) * prob_sv(:, j)') .* reshape(reshape(Kts, ntest * nsv, nkernel) * betas(j, :)', ntest, nsv);
end
s = kaux * model.coefsup + model.w0;
end
