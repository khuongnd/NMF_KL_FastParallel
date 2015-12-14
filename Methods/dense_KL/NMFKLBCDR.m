function [W, F, his] = NMFKLBCDR(V, K, opts),
    if issparse(V),
        V = full(V);
    end
    [W, F, his.errors, his.times] = NMFKLDenseR(V, K, opts.maxIter, opts.W0', opts.F0, opts.verbose);
    W = W';
end