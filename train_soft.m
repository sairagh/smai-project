function [model, history0, status] = train_soft(K, y, prob, C, mkl_norm, options)

options = set_defaults(options,...
                               'nbitermax', 100, ...
                               'epssvm', 1e-2, ...                        
                               'epsbeta', 1e-4, ...                        
                               'epsdualgap', 1e-3, ...                      
                               'verbose_last', 1, ...
                               'verbose', 0);

alphainit = [];                                                            
[,~,nkernel] = size(K);                                               
nclust = size(prob, 2);                                                
nclass = length(unique(y));
if nclass == 1                   
    model.w0 = y(1);
    model.indsup = 1;
    model.coefsup = 0;
    model.prob_sv = zeros(1, nclust);
    model.betas = zeros(nclust, nkernel);
    history0.obj = [];
    status = 1;
    return;
end
if nclass > 2
    error('this is not a two-class problem');
end


betas = ones(nclust, nkernel) / (nkernel^(1/mkl_norm));                    

dualSquare = zeros(nclust, nkernel);                                        

status = 0;                                                             
history0.obj_dual = [];
history0.dual_gap = [];
history0.obj_prime = [];
history0.eps = [];
iter = 0;                                                                  
obj_dual = -inf;
while iter < options.nbitermax && ~status
    iter = iter + 1;
    betasOld = betas;
    obj_dual_old = obj_dual;
X	    Kaux = zeros(ndata, ndata);
    for j = 1 : nclust              
        Kaux = Kaux + (prob(:, j) * prob(:, j)') .* reshape(reshape(K, ndata * ndata, nkernel) * betas(j, :)', ndata, ndata);
    end
    [indsup, coefsup, w0, obj_prime, alphainit] = svm_solve(Kaux, y, C, options, alphainit);
  
    
    for j = 1 : nclust
        tmp = coefsup .* prob(indsup, j);   
        for m = 1 : nkernel            
            dualSquare(j, m) = tmp' * K(indsup, indsup, m) * tmp;            
        end
    end        
    weightSquare = betas.^2 .* dualSquare;    
    
    ty_pred = Kaux(:, indsup) * coefsup + w0;                              
    obj_prime = sum(weightSquare(:)./(2 * betas(:)+eps)) + C * sum(max(1 - ty_pred .* y, 0)); 
    obj_dual = 0;                                                             
    for j = 1 : nclust
        sumWeight = sum(weightSquare(j, :) .^ (mkl_norm/(mkl_norm+1))); 
        sumWeight = max(sumWeight, eps);                                           
        betas(j, :) = (weightSquare(j, :).^(1/(mkl_norm+1))) / (sumWeight^(1/mkl_norm));        
        if mkl_norm == 1
            obj_dual = obj_dual + max(dualSquare(j, :));
        else
            sumdual = sum(dualSquare(j, :) .^ (mkl_norm/(mkl_norm-1))); 
            obj_dual = obj_dual + sumdual ^ ((mkl_norm-1)/mkl_norm);
        end
    end
    obj_dual = sum(coefsup .* y(indsup)) - obj_dual / 2;    
    dual_gap = abs(1 - obj_dual / obj_prime);
    dif_obj = abs(1 - obj_dual / obj_dual_old);    
    dif_beta = norm(betasOld - betas, 'inf');
    
    if dual_gap < options.epsdualgap || dif_beta < options.epsbeta || dif_obj < 1e-6
        status = 1;
    end
    options.epssvm = min([options.epssvm, 1e-2 * dual_gap, 5e-2 * dif_beta, 5e-2 * dif_obj]);            
    options.epssvm = max(options.epssvm, 1e-6);    
    if iter > options.nbitermax * 0.66
        options.epssvm = 1e-6;
    end
    history0.obj_prime = [history0.obj_prime obj_prime];
    history0.obj_dual = [history0.obj_dual obj_dual];
    history0.dual_gap = [history0.dual_gap dual_gap];
    history0.eps = [history0.eps, options.epssvm];
    if options.verbose == 1
        fprintf('iter = %d obj = %5.4f\n', iter, obj_dual);
    end          
end

model.coefsup = coefsup;
model.w0 = w0;
model.indsup = indsup;
model.prob_sv = prob(indsup, :);
model.betas = betasOld;

if ~status && options.verbose_last
	dual_gap = history0.dual_gap(end)
	dif_beta
	epss = options.epssvm(end)
    fprintf('the maximum allowed iterations has been reached (CLMKL)\n');
end

