function opts = getOptions(V, K)
    [N, M] = size(V);
    opts.maxIter = 10;
    opts.tolerance = 1e-3;
    opts.threads = 1;
    opts.verbose = 1;
    opts.maxThread = 1;
    opts.params = [0, 0, 0, 0]';
end