function z=teststat_k4(X0_ord_rep,X1_rep,pol,o ,u,p,test)
%TESTSTAT_K4 returns z=teststatistic of location polynomial at locations 0:1/(p-1):1
% and the deviations from quantile bounds for two dimensional data with location polynomials of degree four. Interpolation between copula samples X0 and X1
% Inputs:
%   X0_ord- ordered sample X0, s.t the minimal mathching is between corresponding X0_ord(i) and X1(i). X0_ord repeated p times
%   X1- original copula sample, repeated p times
%   pol- 3*(original sample size) parameters of location polynomials + penalized deviations from quantile bounds 
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
m=3*n/p;
k=1/(p-1);
   trA=@(x)12./(1+exp(-x))-6;
   trB=@(a,b)(2.*sqrt(2).*sqrt(36-trA(a).^2))./(6*(1+exp(-b)))-(sqrt(2)*sqrt(36-trA(a).^2)+4.*trA(a))./6;
   trC=@(a,b,c)(2*sqrt(36-9.*trA(a).^2-24.*trA(a).*trB(a,b)-18.*trB(a,b).^2))./(1+exp(-c))-sqrt(36-9.*trA(a).^2-24.*trA(a).*trB(a,b)-18.*trB(a,b).^2);

   t4d=@(a,b,c)(sqrt(2)*sqrt(-9*a^2-24*a*b-18*b^2-c^2+36)-2*c)./6;
   t4a=@(a,b,c)t4d(a,b,c).^2;
   t4b=@(a,b,c)(b.^2-t4d(a,b,c).^2+2*c*t4d(a,b,c))./2;
   t4c=@(a,b,c)(2*a*b-2*c*t4d(a,b,c)+c^2)./3;
   p4=@(pol)(ones(m/3,1)-arrayfun(t4a,trA(pol(1:m/3)),trB(pol(1:m/3),pol(m/3+1:2*m/3)),trC(pol(1:m/3),pol(m/3+1:2*m/3),pol(2*m/3+1:3*m/3)))-arrayfun(t4b,trA(pol(1:m/3)),trB(pol(1:m/3),pol(m/3+1:2*m/3)),trC(pol(1:m/3),pol(m/3+1:2*m/3),pol(2*m/3+1:3*m/3)))-arrayfun(t4c,trA(pol(1:m/3)),trB(pol(1:m/3),pol(m/3+1:2*m/3)),trC(pol(1:m/3),pol(m/3+1:2*m/3),pol(2*m/3+1:3*m/3)))).*(0:k:1).^4+arrayfun(t4c,trA(pol(1:m/3)),trB(pol(1:m/3),pol(m/3+1:2*m/3)),trC(pol(1:m/3),pol(m/3+1:2*m/3),pol(2*m/3+1:3*m/3))).*(0:k:1).^3+arrayfun(t4b,trA(pol(1:m/3)),trB(pol(1:m/3),pol(m/3+1:2*m/3)),trC(pol(1:m/3),pol(m/3+1:2*m/3),pol(2*m/3+1:3*m/3))).*(0:k:1).^2+arrayfun(t4a,trA(pol(1:m/3)),trB(pol(1:m/3),pol(m/3+1:2*m/3)),trC(pol(1:m/3),pol(m/3+1:2*m/3),pol(2*m/3+1:3*m/3))).*(0:k:1);
pol2=pol(m+1:m+4*(p-2));
pol=pol(1:m);
as= reshape(p4(pol),n,1);
alphas=[as as];
X2= (ones(n,2)-alphas).*X0_ord_rep+alphas.*X1_rep;

z=[];
for i=1:(p-2) %do not test margins of X0 and X1
   z=[z teststat(X2(m/3*i+1:m/3*i+m/3,1),m/3,test)-o-pol2(1+4*(i-1))];
   z=[z u-teststat(X2(m/3*i+1:m/3*i+m/3,1),m/3,test)-pol2(2+4*(i-1));];
z=[z teststat(X2(m/3*i+1:m/3*i+m/3,2),m/3,test)-o-pol2(3+4*(i-1))];
z=[z u-teststat(X2(m/3*i+1:m/3*i+m/3,2),m/3,test)-pol2(4+4*(i-1))];
end
z=[z -pol2' ];