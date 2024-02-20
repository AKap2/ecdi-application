function T1=teststat(X2,n, test)
%TESTSTAT returns the teststatistic for a given test 
% it can be extended if different teststatistics are used. 
% 'CS-A' is the teststatistic from the corresponding paper. 
% 'T2' - 'T5' are alternative teststatistics for uniformity that are differentiable
% k - scaling factor that minimizes the approximation error of the range estimator approx_max()-approx_min()
% 1/m - expected value of range estimator approx_max()-approx_min() under assumption X2 \sim U(0,1)

%For k =[1,2,5,10,20,50,100,700], 1/m=[13.8779456864436	14.1178565308825	15.5609754326584	19.1852364964591	27.7939382142961	55.9464582331057	104.534656275689	700.249701087256];

switch test
    case 'T2'
        k=10; 
        m=1/19.1852364964591; 
        approx_min=@(Data) -log(sum(exp(-Data)));
        approx_max=@(Data) log(sum(exp(Data)));
        T1=(( 1/k*(approx_max(k*X2)-approx_min(k*X2))*m   )^2  ) / ( sum((X2-0.5).^2) );

    case 'T3'
        T1=mean(X2)/sum((X2-0.5.*ones(n,1)).^2);

    case 'T4'
        k=1; 
        approx_min=@(Data) -log(sum(exp(-Data)));
        approx_max=@(Data) log(sum(exp(Data)));
        T1=(( 1/k*(approx_max(k*X2)-approx_min(k*X2))/mean(X2) )^2  ) / ( sum((X2-0.5).^2) );

    case 'T5'
        k=1; 
        m=1/13.8779456864436; 
        approx_min=@(Data) -log(sum(exp(-Data)));
        approx_max=@(Data) log(sum(exp(Data)));
        T1=(( 1/k*(approx_max(k*X2)-approx_min(k*X2))*m   )^2  ) / ( (n-1)*var(X2)*(2-1/n*sum(((X2-mean(X2))./sqrt(var(X2))).^3))^2 );


    case 'CS-A'
        k=1; 
        m=1/13.8779456864436; 
        approx_min=@(Data) -log(sum(exp(-Data)));
        approx_max=@(Data) log(sum(exp(Data)));
        T1=(( (approx_max(k*X2)-approx_min(k*X2))*m   )^2 ) / ( (n-1)*var(X2) );
end