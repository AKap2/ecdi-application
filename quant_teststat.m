function [z,u,o] = quant_teststat(n,m,proc,test)
%QUANT_TESTSTAT simulates the 'proc'/2 and 1-'proc'/2 quantiles of 'test'.
%Inputs:
%   n - required sample size
%   m - repetitions of sampling
%   proc -  alpha level of the two-tailed test  
%   test - used teststatistic ('CS-A', 'T2'-'T5'implemented, for others alter teststat.m) 
%Outputs:
%   z - m teststatistics 'test' for uniform samples of size n, ordered ascendingly
%   u - 'proc'/2 quantile of 'test' for samples of size 'n'
%   o - 1-'proc'/2 quantile of 'test' for samples of size 'n' 
%
%Author: Anika Kaplan
z=ones(m,1);
rand('seed',123);
W=rand(n,m);
for i=1:m
    z(i)= teststat(W(:,i), n,test);
end
z=sort(z);
u=z(round(m*proc/2));
o=z(round(m*(1-proc/2)));