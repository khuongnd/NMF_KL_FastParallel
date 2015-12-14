methods = {@NMFKLExact, @nmf_kl, @NMFKLBCD};
%methods = {@nmf_kl};
names = {'SRCCD', 'MU', 'CCC'};
datasets = {'digits', 'Reuters21578', 'TDT2', 'RCV1_4Class'};
%datasets = {'digits', 'Reuters21578', 'TDT2'};


set(groot,'defaultAxesLineStyleOrder',{'-*',':','o'}) 
styles = {'r-o', 'g-d', 'm:*', 'y-.+'};
plts = [];
K = 10;
for d=1:length(datasets),
    subplot(1, length(datasets), d)
    hold on;
%     minError = 1e100;
%     for i=1:length(methods),
%         load(sprintf('result_all/%s_%d_%s.mat', datasets{d}, K, func2str(methods{i})));
%         minError = min(minError, min(errors));
%     end
    for i=1:length(methods),
        load(sprintf('result_all/%s_%d_%s.mat', datasets{d}, K, func2str(methods{i})));
        loglog(times, errors, styles{i}, 'MarkerSize',3);
        %folder = '~/Dropbox/CRCDForNMFKL/figs';
        folder = '~/Dropbox/PhDThesis/chapters/05_figs';
        datafile = sprintf('%s/%s_%d_%s.data', folder, datasets{d}, K, func2str(methods{i}));
        fprintf('%f %s\n', errors(length(errors)), datafile);
        %size(times)
        %errors = (errors - minError) / minError;
        if (size(times, 1) == 1),
            dlmwrite(datafile, [times'  errors'], ' ');
        else
            dlmwrite(datafile, [times  errors], ' ');
        end;
    end
    %plts = [plts get(gca,'Children')];
    legend(plts, names{1}, names{2}, names{3}, 'Location','northoutside','Orientation','horizontal');
end
shg();   
