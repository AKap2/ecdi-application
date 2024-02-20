[Xt,Test1_2,Test1_3,Test1_4,Test2_2,Test2_3,Test2_4, alpha, A,figLocationPoly,figTest2, figTest3,figTest4]=ecdi2dparam(family0,family1,alpha0,alpha1,n,test,iterations,mu,quant,t,p,KK)
%ECDI2DPARAM creates copula sample X(t) at t \in (0,1) between X(0) and X(1), evaluations and auxiliary results
%Input: 
%   family0- type of parametric copula 1 ('Gauss','Clayton','Gumbel','Frank','t',...)
%   family1 - type of parametric copula 2
%   alpha0 - parameter of copula 1 (corresponding to copula type)
%   alpha1 - parameter of copula 2 (corresponding to copula type)
%   n - size of copula samples
%   test - teststatistic for uniformity (differentiable)
%           ('CS-A','T2','T3','T4','T5')    
%   iterations - \neq 1, if several samples should be produced from different initial copula samples with the same parameters
%   mu - scales penalty term (=10000)
%   quant - quantifies the level of the used teststatistic (\in (0,0.5))
%   t - position of the interpoated sample (t \in [0,1])
%   p - Number of evaluation points including t=0 and t=1
%   KK - maximum degree of location polynomials
%Output: 
%   Xt - copula sample at interpolation position t in [0,1] between X0 and X1 from given copula families
%   Test1_2 - Value of teststatistic 'test' of samples at positions [0:1/(p-1):1], first dimension, KK=2 
%   Test1_3 - Value of teststatistic 'test' of samples at positions [0:1/(p-1):1], first dimension, KK=3
%   Test1_4 - Value of teststatistic 'test' of samples at positions [0:1/(p-1):1], first dimension, KK=4
%   Test2_2 - Value of teststatistic 'test' of samples at positions [0:1/(p-1):1], second dimension, KK=2
%   Test2_3 - Value of teststatistic 'test' of samples at positions [0:1/(p-1):1], second dimension, KK=3
%   Test2_4 - Value of teststatistic 'test' of samples at positions [0:1/(p-1):1], second dimension, KK=4
%   alpha - Parameters of location polynomial
%   A - values of location polynomial at positions  [0:1/(p-1):1]
%   figLocationPoly - plot of location polynomial (spline based on A)
%   figTest2 - Values Test1_2 and Test2_2 displayed in boxplot
%   figTest3 - Values Test1_3 and Test2_3 displayed in boxplot
%   figTest4 - Values Test1_4 and Test2_4 displayed in boxplot
%
% Author: Anika Kaplan


k=1/(p-1);% Interpolation step size, RESTRICTION - to implement: 1:k being an integer
P=[0:k:1]; %Interpolation parameter    
[~,u,o]=quant_teststat(n,mu,quant,test); %Upper and Lower bound for U(0,1) test statistic, Sample size n, a/2 und 1-a/2 quantiles
 
  %Store results
    Test1_2=[];
    Test1_3=[];
    Test1_4=[];
    Test2_2=[];
    Test2_3=[];
    Test2_4=[];
    w_save=[];
    for iters = 1:iterations
      rng(iters+100,'twister')
      %Create samples between which we interpolate   
        X0=copularnd(family0, alpha0, n);
        X1=copularnd(family1, alpha1, n);
        %generate copula samples
       while (teststat(X0(:,1),n,test)>o)||(teststat(X0(:,1),n,test)<u)||(teststat(X0(:,2),n,test)>o)||(teststat(X0(:,2),n,test)<u)
            X0=copularnd(family0, alpha0, n);
        end
        while (teststat(X1(:,1),n,test)>o)||(teststat(X1(:,1),n,test)<u)||(teststat(X1(:,2),n,test)>o)||(teststat(X1(:,2),n,test)<u)
            X1=copularnd(family1, alpha1, n);
        end
      

        %% Linear Sum Assignment Problem via Hungarian Method , minimum matching X0 and X1
          D=distancematrix(X0,X1);
        [V1]=matchpairs(D,1000,'min');

    %% Apply Matching
        X0_ord=zeros(n,2);
        w=zeros(n,1);
        for i=1:n
         X0_ord(i,:)=X0(V1(i),:);
         w(i)=norm(X1(i,:)-X0_ord(i,:));
        end
        %Length of the matching edges (=distance between X0 and X1)
        W=sum(w);
        w_save=[w w_save];
               
    %% Minimization via fmincon.margins are U(0,1)
    p=length(P);
            w=repmat(w,p,1);
            X0_ord_rep=repmat(X0_ord,p,1);
            X1_rep=repmat(X1,p,1);
            P_neu=repelem((0:k:1)',1,n);
    if KK==2
        logis = @(x)2./(1+exp(-x));
        x0=zeros(n+4*(p-2),1); % Start with b=0.

        options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'ConstraintTolerance', 10^(-6), 'StepTolerance', 10^(-25));
        problem.options = options;
        problem.solver = 'fmincon';
        problem.objective= @(pol) sum((sum((reshape(w,n,p).*((0:k:1).*logis(pol(1:n))+(0:k:1).^2.*(ones(n,1)-logis(pol(1:n)))-P)))).^2)+ mu*sum(pol(n+1:n+4*(p-2)));
        problem.x0 = x0;
        problem.nonlcon = @(pol)teststat_k2(X0_ord_rep ,X1_rep ,pol ,o ,u, p,test);
        [alpha,f,exitflag] = fmincon(problem);
        alpha(1:n)=logis(alpha(1:n));
        alpha=alpha(1:n);
        A= reshape(((0:k:1).*(alpha(1:n))+(0:k:1).^2.*(ones(n,1)-alpha(1:n))),n*length((0:k:1)),1);%alphas finden
        alphas=[A A];
        X2= (ones(n*p,2)-alphas).*X0_ord_rep+alphas.*X1_rep;
        %
        %Save Test statistics
            d=0;
            for j=0:k:1 
                Test1_2=[Test1_2 teststat(X2(n*d+1:n*d+n,1),n,test)];
                d=d+1;
            end
            d=0;
            for j=0:k:1 
                Test2_2=[Test2_2 teststat(X2(n*d+1:n*(d+1),2),n,test)];
                d=d+1;
            end
        A_t=t.*alpha(1:n)+t^2.*(ones(n,1)-alpha(1:n));
        alpha_t=[A_t A_t];
        Xt=(ones(n,2)-alpha_t).*X0_ord+alpha_t.*X1;
    end

    if KK==3
        x0=[zeros(n,1); 100000.*ones(n,1); zeros(4*(p-2),1)];
        trA3=@(x)2.*sqrt(12)./(1+exp(-x))-sqrt(12);
        trB3=@(a,b)(sqrt(3).*sqrt(12-trA3(a).^2)-3.*trA3(a))./(3*(1+exp(-b)))-(sqrt(3)*sqrt(12-trA3(a).^2)-3.*trA3(a))./6;
        t3a=@(b) b^2;
        t3c=@(a,b) -2*a^2-6*a*b-6*b^2+6;
        t3b=@(a,b)a*b+t3c(a,b)/2;
        p3=@(pol) (ones(n,1)-arrayfun(t3a,trB3(pol(1:n),pol(n+1:2*n)))-arrayfun(t3b,trA3(pol(1:n)),trB3(pol(1:n),pol(n+1:2*n)))).*(0:k:1).^3+arrayfun(t3b,trA3(pol(1:n)),trB3(pol(1:n),pol(n+1:2*n))).*(0:k:1).^2+arrayfun(t3a,trB3(pol(1:n),pol(n+1:2*n))).*(0:k:1);
        options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'ConstraintTolerance', 10^(-6), 'StepTolerance', 10^(-5));
        problem.options = options;
        problem.solver = 'fmincon';
        problem.objective= @(pol) sum((sum((reshape(w,n,p).*(p3(pol) -P)))).^2)+ mu*sum(pol(2*n+1:2*n+4*(p-2)));
        problem.x0 = x0;
        problem.nonlcon = @(pol)teststat_k3(X0_ord_rep ,X1_rep ,pol ,o ,u, p,test);
        [alpha,f,exitflag] = fmincon(problem);
        alpha=alpha(1:2*n);
        A= reshape(p3(alpha),n*p,1);
        alphas=[A A];
        X2= (ones(n*p,2)-alphas).*X0_ord_rep+alphas.*X1_rep;
        %Save Test statistics
            d=0;
            for j=0:k:1 
                Test1_3=[Test1_3 teststat(X2(n*d+1:n*(d+1),1),n,test)];
                d=d+1;
            end
            d=0;
            for j=0:k:1 
                Test2_3=[Test2_3 teststat(X2(n*d+1:n*(d+1),2),n,test)];
                d=d+1;
            end
        A_t=(ones(n,1)-arrayfun(t3a,trB3(alpha(1:n),alpha(n+1:2*n)))-arrayfun(t3b,trA3(alpha(1:n)),trB3(alpha(1:n),alpha(n+1:2*n)))).*t^3+arrayfun(t3b,trA3(alpha(1:n)),trB3(alpha(1:n),alpha(n+1:2*n))).*t^2+arrayfun(t3a,trB3(alpha(1:n),alpha(n+1:2*n))).*t;
        alpha_t=[A_t A_t];
        Xt=(ones(n,2)-alpha_t).*X0_ord+alpha_t.*X1;
    end


    if KK==4
        x0=[zeros(n,1); log(3+2*sqrt(2))*ones(n,1); zeros(n,1); zeros(4*(p-2),1)];
        trA=@(x)12./(1+exp(-x))-6;
        trB=@(a,b)(2.*sqrt(2).*sqrt(36-trA(a).^2))./(6*(1+exp(-b)))-(sqrt(2)*sqrt(36-trA(a).^2)+4.*trA(a))./6;
        trC=@(a,b,c)(2*sqrt(36-9.*(trA(a).^2)-24.*trA(a).*trB(a,b)-18.*(trB(a,b).^2)))./(1+exp(-c))-sqrt(36-9.*trA(a).^2-24.*trA(a).*trB(a,b)-18.*trB(a,b).^2);

        t4d=@(a,b,c)(sqrt(2)*sqrt(-9*a^2-24*a*b-18*b^2-c^2+36)-2*c)./6;
        t4a=@(a,b,c)t4d(a,b,c).^2;
        t4b=@(a,b,c)(b.^2-t4d(a,b,c).^2+2*c*t4d(a,b,c))./2;
        t4c=@(a,b,c)(2*a*b-2*c*t4d(a,b,c)+c^2)./3;
        p4=@(pol)(ones(n,1)-arrayfun(t4a,trA(pol(1:n)),trB(pol(1:n),pol(n+1:2*n)),trC(pol(1:n),pol(n+1:2*n),pol(2*n+1:3*n)))-arrayfun(t4b,trA(pol(1:n)),trB(pol(1:n),pol(n+1:2*n)),trC(pol(1:n),pol(n+1:2*n),pol(2*n+1:3*n)))-arrayfun(t4c,trA(pol(1:n)),trB(pol(1:n),pol(n+1:2*n)),trC(pol(1:n),pol(n+1:2*n),pol(2*n+1:3*n)))).*(0:k:1).^4+arrayfun(t4c,trA(pol(1:n)),trB(pol(1:n),pol(n+1:2*n)),trC(pol(1:n),pol(n+1:2*n),pol(2*n+1:3*n))).*(0:k:1).^3+arrayfun(t4b,trA(pol(1:n)),trB(pol(1:n),pol(n+1:2*n)),trC(pol(1:n),pol(n+1:2*n),pol(2*n+1:3*n))).*(0:k:1).^2+arrayfun(t4a,trA(pol(1:n)),trB(pol(1:n),pol(n+1:2*n)),trC(pol(1:n),pol(n+1:2*n),pol(2*n+1:3*n))).*(0:k:1);
        options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'ConstraintTolerance', 10^(-6), 'StepTolerance', 10^(-6));
        problem.options = options;
        problem.solver = 'fmincon';
        problem.objective= @(pol) sum((sum((reshape(w,n,p).*(p4(pol) -P)))).^2)+ mu*sum(pol(3 *n+1:3*n+4*(p-2)));
        problem.x0 = x0;
        problem.nonlcon = @(pol)teststat_k4(X0_ord_rep ,X1_rep ,pol ,o ,u, p,test);
        [alpha,f,exitflag] = fmincon(problem);
        alpha=alpha(1:3*n);
        A= reshape(p4(alpha(1:3*n)),n*p,1);
        alphas=[A A];
        X2= (ones(n*p,2)-alphas).*X0_ord_rep+alphas.*X1_rep;

    %Accumulate Test statistics (needed when several iterations)
        d=0;
        for j=0:k:1 
            Test1_4=[Test1_4  teststat(X2(n*d+1:n*(d+1),1),n,test)];
            d=d+1;
        end
        d=0;
        for j=0:k:1 
            Test2_4=[Test2_4 teststat(X2(n*d+1:n*(d+1),2),n,test)];
            d=d+1;
        end
        A_t=(ones(n,1)-arrayfun(t4a,trA(alpha(1:n)),trB(alpha(1:n),alpha(n+1:2*n)),trC(alpha(1:n),alpha(n+1:2*n),alpha(2*n+1:3*n)))-arrayfun(t4b,trA(alpha(1:n)),trB(alpha(1:n),alpha(n+1:2*n)),trC(alpha(1:n),alpha(n+1:2*n),alpha(2*n+1:3*n)))-arrayfun(t4c,trA(alpha(1:n)),trB(alpha(1:n),alpha(n+1:2*n)),trC(alpha(1:n),alpha(n+1:2*n),alpha(2*n+1:3*n)))).*t^4+arrayfun(t4c,trA(alpha(1:n)),trB(alpha(1:n),alpha(n+1:2*n)),trC(alpha(1:n),alpha(n+1:2*n),alpha(2*n+1:3*n))).*t^3+arrayfun(t4b,trA(alpha(1:n)),trB(alpha(1:n),alpha(n+1:2*n)),trC(alpha(1:n),alpha(n+1:2*n),alpha(2*n+1:3*n))).*t^2+arrayfun(t4a,trA(alpha(1:n)),trB(alpha(1:n),alpha(n+1:2*n)),trC(alpha(1:n),alpha(n+1:2*n),alpha(2*n+1:3*n))).*t;
       alpha_t=[A_t A_t];
        Xt=(ones(n,2)-alpha_t).*X0_ord+alpha_t.*X1;  
end



%Values at p evaluated points 
A=reshape(A(1:n*length(P)),[n, length(P)]);


    %% Spline Interpolation of the alphai s, boxplots
    if iters ==1
        figLocationPoly =figure(1);
        x = [0:k:1];
        xq1 = 0:0.1*k:1;
        for kk=1:n
            y = [ A(kk,:) ]; 
            pch = pchip(x,y,xq1);
            plot(x,y,'o',xq1,pch,'-'); 
            ylim([0 1])
            hold on; 
        end
        hold off; 
    end;
end;

   
figTest2=figure;
    W=[Test1_2;Test2_2];
    line([1,21],[u,u]);hold on;
    line([1,21],[o,o]);hold on;
    boxplot(W);
    xlim([0,p+1])
    xlabel('Location in [0:k/21:1]')
    ylabel('Value of test statistic')

figTest3=figure;
    W=[Test1_3;Test2_3];
    line([1,21],[u,u]);hold on;
    line([1,21],[o,o]);hold on;
    boxplot(W);
    xlim([0,p+1])
    xlabel('Location in [0:k/21:1]')
    ylabel('Value of test statistic')
        
figTest4=figure;
    W=[Test1_4;Test2_4];
    line([1,21],[u,u]);hold on;
    line([1,21],[o,o]);hold on;
    boxplot(W);
    xlim([0,p+1])
    xlabel('Location in [0:k/21:1]')
    ylabel('Value of test statistic')
end