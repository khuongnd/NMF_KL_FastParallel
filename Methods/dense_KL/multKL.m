function [W,H,objlist,timelist] = multKL(V,k,maxiter,W0,H0)

n = size(V,1);
m = size(V,2);

W= W0(:,1:k);
H = H0(1:k,:);
total = 0;

for iter = 1:maxiter
	begin=cputime;
	W = W.*( (((V+1e-5)./(W*H+1e-5))*H') ./repmat(sum(H,2)',n,1));
	H = H.*( (W'*((V+1e-5)./(W*H+1e-5))) ./repmat(sum(W,1)',1,m));
	total = total + cputime-begin;
	timelist(iter) = total;
	objlist(iter) = KL(V,W,H);
	fprintf('Iter %g obj %g\n', iter, objlist(iter));
end
