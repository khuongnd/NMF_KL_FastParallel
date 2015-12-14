addpath (genpath('Methods'));
graphics = 0;
if (graphics),
    hold on;
end;

maxThreads = [1];

methods = {@nmf_kl_sparse_v, @nmf_kl, @NMFKLBCD, @NMFKLExact, @NMFKLAppro};
his = {[], [], []};

datasets = {'RCV1_4Class'};

methods = {@NMFKLExactPar};

set(groot,'defaultAxesLineStyleOrder',{'-*',':','o'}) 
styles = {'r-o', 'g-d', 'k:*', 'y-.+'};

for maxThread=maxThreads,
    for i=1:length(methods),
        for d=1:length(datasets),
            K = 10;
            fprintf('%s, %s, K=%d\n', datasets{d}, func2str(methods{i}), K);
            load(sprintf('datasets/%s_orgin.mat', datasets{d}));
            load(sprintf('datasets/%s_%d.mat', datasets{d}, K));
            opts = getOptions(V, K);
            opts.maxIter = 10;
            opts.W0 = W0;
            opts.F0 = F0;
            opts.maxThread = maxThread;
            %opts.params = [0, 0.01, 0, 0.01]';
            %opts.params = [0, 0.0, 0, 0.0]';
            [W, F, his{i}] = methods{i}(V, K, opts);
            if (graphics),
                loglog(his{i}.times, his{i}.errors, styles{i}, 'MarkerSize',3);
            end;
            times = his{i}.times;
            errors = his{i}.errors;
            wSparse = 1.0 - 1.0*sum(sum(W~=0))/sum(sum(W==W))
            fSparse = 1.0 - 1.0*sum(sum(F~=0))/sum(sum(F==F))
            %outFile = sprintf('results/%s_%d_%s_%d.mat', datasets{d}, K, func2str(methods{i}), maxThread);
            %save(outFile, 'W', 'F', 'times', 'errors', 'wSparse', 'fSparse');
        end
    end
end;

