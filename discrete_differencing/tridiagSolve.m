function [y]=tridiagSolve(a,b,c, d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: thomas
%
% To solve a tridiagonal
% linear system using the
% Thomas algorithm
%
% The ith equation is:
%
% a(i)*y(i-1) + b(i)*y(i) + c(i)*y(i+1) = d(i) 
%
% N x N matrix
%
% note that a(1)=0 and c(N)=0
%
% Andrew Duggleby 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(d);

for k=2:N
    m=a(k)/b(k-1);
    b(k)=b(k)-m*c(k-1);
    d(k)=d(k)-m*d(k-1);
end
y(N)=d(N)/b(N);

for k=N-1:-1:1
    y(k)=(d(k)-c(k)*y(k+1))/b(k);
end
