clear all;close all;
P = [1, 1000];
Lambda = [1e-4, 1e-6];
m = 10;
sigma = 40;
s_lim1 = 1;
s_lim2 = 10;
n_div = 10;
x_arr = [];
s_arr = [];

%parfor s=linspace(s_lim1,s_lim2,n_div)
%    x = ComputeRateThomasType1user(P,Lambda,m,sigma,s)
%    x_arr = [x_arr x];
%    s_arr = [s_arr s];
%end

s=10000:1000:20000;
parfor i=1:length(s)
 x= ComputeRateThomasType1user(P,Lambda,m,sigma,s(i))
 x_arr = [x_arr x];
 s_arr = [s_arr s(i)];
end
s_arr
x_arr
