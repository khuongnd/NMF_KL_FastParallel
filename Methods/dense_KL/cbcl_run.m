clear
max_iter = [100 1200];
rand('seed',0);
load '../data/cbcl';
n=size(V,1);
m=size(V,2);
k=49;

times = 1;
objlist0 = zeros(1,max_iter(1));
objlist1 = zeros(1,max_iter(2));
timelist0 = zeros(1,max_iter(1));
timelist1 = zeros(1,max_iter(2));
for i=1:times
	W = abs(rand(n,k));
	H = abs(rand(k,m));
	obj_ini = KL(V,W,H);
	[w h obj0 time0] = KLnmf(V,k,max_iter(1),W',H, 1);
	[w1 h1 obj1 time1] = multKL(V,k,max_iter(2), W, H );
	objlist0 = objlist0 + obj0;
	objlist1 = objlist1 + obj1;
	timelist0 = timelist0 + time0;
	timelist1 = timelist1 + time1;
end
objlist0 = objlist0/times;
objlist1 = objlist1/times;
timelist0 = timelist0/times;
timelist1 = timelist1/times;

