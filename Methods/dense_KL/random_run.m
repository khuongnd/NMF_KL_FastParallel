rand('seed',0);
max_iter = [100 1200];

%%% statistics of synthetic data
m = 500;
n = 1000;
k = 10;
rate = 0.3;

%%% generate synthetic data
W_org = rand(n,k);, W_org(rand(n,k)<rate)=0;
H_org = rand(k,m);, H_org(rand(k,m)<rate)=0;
V = W_org * H_org;

%%% number of random initial points
times = 1

%%% initialize 
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
