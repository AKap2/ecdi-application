function z=teststat_k2_3d(X0_ord,X1,pol,o ,u,p,test)
%TESTSTAT_K2_3d returns z=teststatistic of location polynomial at locations 0:1/(p-1):1
% and the deviations from quantile bounds for three dimensional data with location polynomials of degree two. Interpolation between copula samples X0 and X1
% Inputs:
%   X0_ord- ordered sample X0, s.t the minimal mathching is between corresponding X0_ord(i) and X1(i). X0_ord repeated p times
%   X1- original copula sample, repeated p times
%   pol- (original sample size) parameters of location polynomials + penalized deviations from quantile bounds 
%   o- upper quantile bound of used teststatistic
%   u- lower quantile bound of used teststatistic
%   p- number of points, where the location polynomial is evaluated including X0 and X1
% 
%Output:
%   z- teststatistics for all three margins of each evaluation point (1/(p-1) to (p-2)/(p-1))
%      and the corresponding deviations (if above o or below u)
%
%Author: Anika Kaplan

n=length(X0_ord);
m=n/p;
logis = @(x)2./(1+exp(-x));
pol2=pol(m+1:m+6*(p-2),1);
pol=logis(pol(1:m));

as= reshape(((0:1/(p-1):1).*pol+(0:1/(p-1):1).^2.*(ones(m,1)-pol)),n,1);
alphas=[as as as];
X2= (ones(n,3)-alphas).*X0_ord+alphas.*X1;


z=[];
for i=1:(p-2) %do not test margins of X0 and X1
    z=[z teststat(X2(m*i+1:m*i+m,1),m,test)-o-pol2(1+6*(i-1))];
   z=[z u-teststat(X2(m*i+1:m*i+m,1),m,test)-pol2(2+6*(i-1));];
z=[z teststat(X2(m*i+1:m*i+m,2),m,test)-o-pol2(3+6*(i-1))];
z=[z u-teststat(X2(m*i+1:m*i+m,2),m,test)-pol2(4+6*(i-1))];
z=[z teststat(X2(m*i+1:m*i+m,3),m,test)-o-pol2(5+6*(i-1))];
z=[z u-teststat(X2(m*i+1:m*i+m,3),m,test)-pol2(6+6*(i-1))];
end
z=[z -pol2'];
