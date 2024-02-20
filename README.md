# ecdi-application
 The ECDI method implemented for 2- and 3- dimensional data

 ECDI- A non-parametric copula interpolation with applications in a combustion process and climate data

==================================================================

 ECDI is the empirical copula displacement interpolation. Based on two samples (with uniform margins) from the same environment with only one different external factor $t$ (scaled to [0,1]), additional samples from the same environment are generated with $t \in (0,1)$.

Main functions are:    
ECDI2dparametric - implemented for 2-dimensional copula samples from parametric copulas.     
ECDI3d - implemented for 3-dimensional samples.     

Auxiliary functions are:
distancematrix(X,Y) - function returns a matrix of pairwise distances between elements of samlpe X and sample Y
quant_teststat(n,m,proc,test) - simulates the 'proc'/2 and 1-'proc'/2 quantiles of 'test'
teststat(X2,n, test) - returns the teststatistic of X2 for a given test
teststat_k2(X0_ord,X1,pol,o ,u,p,test) - returns z=teststatistic of location polynomial of degree 2 at locations 0:1/(p-1):1 and deviation from quantiles for 2d data
teststat_k2_3d(X0_ord,X1,pol,o ,u,p,test) - returns z=teststatistic of location polynomial of degree 2 at locations 0:1/(p-1):1 and deviation from quantiles for 3d data
teststat_k3(X0_ord,X1,pol,o ,u,p,test) - returns z=teststatistic of location polynomial of degree 3 at locations 0:1/(p-1):1 and deviation from quantiles for 2d data
teststat_k3_3d(X0_ord,X1,pol,o ,u,p,test) - returns z=teststatistic of location polynomial of degree 3 at locations 0:1/(p-1):1 and deviation from quantiles for 3d data
teststat_k4(X0_ord,X1,pol,o ,u,p,test) - returns z=teststatistic of location polynomial of degree 4 at locations 0:1/(p-1):1 and deviation from quantiles for 2d data
teststat_k4_3d(X0_ord,X1,pol,o ,u,p,test) - returns z=teststatistic of location polynomial of degree 4 at locations 0:1/(p-1):1 and deviation from quantiles for 3d data\\
umarg(X0,k,u,o) - transforms the margins of X0 to be uniform in such a way, that the teststatistic of the margins is larger than u and smaller than o 
