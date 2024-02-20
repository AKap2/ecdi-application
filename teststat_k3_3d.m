function z=teststat_k3_3d(X0_ord_rep,X1_rep,pol,o ,u,p,test)

%TESTSTAT_K3_3d returns z=teststatistic of location polynomial at locations 0:1/(p-1):1
% and the deviations from quantile bounds for three dimensional data with location polynomials of degree three. Interpolation between copula samples X0 and X1
% Inputs:
%   X0_ord- ordered sample X0, s.t the minimal mathching is between corresponding X0_ord(i) and X1(i). X0_ord repeated p times
%   X1- original copula sample, repeated p times
%   pol- 2*(original sample size) parameters of location polynomials + penalized deviations from quantile bounds 
%   o- upper quantile bound of used teststatistic
%   u- lower quantile bound of used teststatistic
%   p- number of points, where the location polynomial is evaluated including X0 and X1
% 
%Output:
%   z- teststatistics for both margins of each evaluation point (1/(p-1) to (p-2)/(p-1))
%      and the corresponding deviations (if above o or below u)
%
%Author: Anika Kaplan

n=length(X0_ord_rep);
m=2*n/p;
k=1/(p-1);
 
   trA3=@(x)2.*sqrt(12)./(1+exp(-x))-sqrt(12);
   trB3=@(a,b)(sqrt(3).*sqrt(12-trA3(a).^2)-3.*trA3(a))./(3*(1+exp(-b)))-(sqrt(3)*sqrt(12-trA3(a).^2)-3.*trA3(a))./6;
   t3a=@(b) b.^2;
   t3c=@(a,b) -2.*a^2-6.*a.*b-6.*b.^2+6;
   t3b=@(a,b)a.*b+t3c(a,b)/.2;
   p3=@(pol) (ones(m/2,1)-arrayfun(t3a,trB3(pol(1:m/2),pol(m/2+1:2*m/2)))-arrayfun(t3b,trA3(pol(1:m/2)),trB3(pol(1:m/2),pol(m/2+1:2*m/2)))).*(0:k:1).^3+arrayfun(t3b,trA3(pol(1:m/2)),(trB3(pol(1:m/2),pol(m/2+1:2*m/2)))).*(0:k:1).^2+arrayfun(t3a,trB3(pol(1:m/2),pol(m/2+1:2*m/2))).*(0:k:1);

pol2=pol(m+1:m+6*(p-2));
pol=pol(1:m);
as= reshape(p3(pol),n,1);
alphas=[as as as];
X2= (ones(n,3)-alphas).*X0_ord_rep+alphas.*X1_rep;

z=[];
for i=1:(p-2) %do not test margins of X0 and X1
   z=[z teststat(X2(m/2*i+1:m/2*i+m/2,1),m/2,test)-o-pol2(1+6*(i-1))];
   z=[z u-teststat(X2(m/2*i+1:m/2*i+m/2,1),m/2,test)-pol2(2+6*(i-1));];
z=[z teststat(X2(m/2*i+1:m/2*i+m/2,2),m/2,test)-o-pol2(3+6*(i-1))];
z=[z u-teststat(X2(m/2*i+1:m/2*i+m/2,2),m/2,test)-pol2(4+6*(i-1))];
z=[z teststat(X2(m/2*i+1:m/2*i+m/2,3),m/2,test)-o-pol2(5+6*(i-1))];
z=[z u-teststat(X2(m/2*i+1:m/2*i+m/2,3),m/2,test)-pol2(6+6*(i-1))];
end
z=[z -pol2']; 