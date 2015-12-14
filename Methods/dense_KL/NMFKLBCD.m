function [W, F, his] = NMFKLBCD(V, K, opts),
    if issparse(V),
        V = full(V);
    end
    [W, F, his.errors, his.times] = NMFKLDense(V, K, opts.maxIter, opts.W0', opts.F0, opts.verbose);
    W = W';
end