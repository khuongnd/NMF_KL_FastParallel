function [obj] = KL(V,W,H)

A = W*H;
obj = sum(sum(V.*log(V+1e-10)-V.*log(A+1e-10)-V+A));
