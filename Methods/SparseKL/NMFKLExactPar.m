function [W, F, his] = NMFKLExactPar(V, K, opts),
    if ~issparse(V),
        V = sparse(V);
    end
    [F, W, his.errors, his.times, his.iters] = NMFSparseKLExactPar(V, opts.F0', opts.W0, opts.maxIter, opts.tolerance, opts.maxThread, opts.params);
    F = F';
end