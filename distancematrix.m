function dm = distancematrix(X,Y)
%DISTANCEMATRIX 
%function returns a matrix of pairwise distances between elements of samlpe X and sample Y
% input sample X=(X1; X2;...) and sample Y=(Y1; Y2,...)
% output matrix dm with dm(i,j)=dist(Xi,Yj)
%Autjor: Anika Kaplan
m= length(X);
n= length(Y);

dm=zeros(m,n);
for i=1:m
    for j=1:n
        dm(i,j)= norm(X(i,:)-Y(j,:));
    end;
end;