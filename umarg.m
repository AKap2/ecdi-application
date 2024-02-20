function X1=umarg(X0,k,u,o)
%UMARG transforms the margins of X0 to be uniform in such a way, that the teststatistic of the margins is larger than u and smaller than o 
% this function is used instead of just ranking the data, to allow more variability
% k is used a a seed 
% Inputs:
%  X0-  of size n and dimension d
%  k - seed for sampling of uniform distribution
%  u - lower bound, chosen alpha/2 quantile of test for samples of size n
%  o - upper bound, chosen 1-alpha/2 quantile of test for samples of size n
% Output: 
%  X1- X0 with transformed margins (overall copula is kept, but positions are altered) 
%Author: Anika Kaplan


n = size(X0,1); % sample size
d = size(X0,2); % dimension
X1 = zeros(n,d);
V1=zeros(n,1);
test='CS-A';
U1=zeros(n,1);
rand('seed',123+k);

for i=1:d
   while (teststat(U1,n,test)>o)||(teststat(U1,n,test)<u)
         U1=rand(n,1);
   end
   [B,rx]=sort(X0(:,i));%/(n+1);
   ru=tiedrank(U1(:,1));%/(n+1);
   V1(ru,:)=U1; %sorted V1 with values from U(0,1)
   X1(rx,i)=V1(:,:);
   U1=zeros(n,1);
end    