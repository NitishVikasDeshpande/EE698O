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

s1=0.001:0.01:0.1;
s2=0.1:0.1:1;
s3=2:2:50;
s4=51:10:100;
s5=100:50:1000;
s6=1001:500:80000;
s=[s1 s2 s3 s4 s5 s6];
fprintf("length of s =%f\n",length(s));
parfor i=1:length(s)
 x= ComputeRateThomasType1user(P,Lambda,m,sigma,s(i))
 x_arr = [x_arr x];
 s_arr = [s_arr s(i)];
end
s_arr
x_arr
Q= trapz(s_arr,x_arr)
fprintf("\n %f",Q)
