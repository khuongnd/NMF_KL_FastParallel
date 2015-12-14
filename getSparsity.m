function res = getSparsity(V)
    res = 1.0*sum(sum(abs(V) == 0)) / sum(sum(V == V)) * 100;
end