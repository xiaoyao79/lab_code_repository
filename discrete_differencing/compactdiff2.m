function [N]=compactdiff2(N,dx,dir)

% if nargin<3
%     %from Lele's original work
     alpha = 2/11; beta=0;    a=12/11;   b=3/11;   c=0;         %tridiagonal,5pt stencil,6th order
%     alpha = 9/38; beta=0;    a=147/152; b=51/95;  c=-23/760;   %tridiagonal,7pt stencil,8th order
%     alpha = 334/899; beta=43/1798; a=1065/1798; b=1038/899; c=79/1798;   %10th order
% end


%BC - Tridiagonal
%i=0, 3rd order
a00 = 13;  a01 = -27;
a02 =  15;  a03 = -1;
alpha01 = 11;
%i=1, 4th order tridiagonal,3pt stencil
alpha10 = 1/10; alpha12 = 1/10;
a10=6/5; a11=-2*6/5; 
a12=6/5; a13 = 0; 
a14 = 0;
%i=2, 6th order
a20 = 3/44; a21 = 12/11;
a22 = -51/22; a23 =  12/11;
a24 =  3/44; a25 = 0;
alpha21 = 2/11; alpha23 = 2/11;

size_N = size(N);

if isempty(dir)
    if size_N(1)~=1
        dir = 1;
    elseif size_N(2)~=1
        dir = 2;
    elseif size_N(3)~=1
        dir = 3;
    else 
        dir = 1;
    end
end

ndim_N = ndims(N);
if ndim_N > 3
    error('function only defined up to 3D matrices')
end

if isempty(dx)
    dx = 1;
end


N = permute(N, [dir, 1:dir-1, dir+1:ndim_N]);
size_N = size(N);   %find again for permuted matrix
ndim_N = ndims(N);


NI = size_N(1);
NJ = size_N(2);
if ndim_N == 3
    NK = size_N(3);
else
    NK = 1;
end
    

if NI < 7
    error('stencil size is 7 pts in active direction')
end

X = zeros(NI,1);
B = zeros(NI,1);
A = sparse([],[],[],NI,NI,3*NI);

RHS = [ c/9 b/4 a -2*(c/9+b/4+a/1) a b/4 c/9 ]/dx.^2;
%LHS = [ beta alpha 1 alpha beta ]; %pentadiagonal
LHS = [ alpha 1 alpha ];  %tridiagonal

%precalculate A matrix 1 time only
        %i = 1;
            A(1,1:2) = [1, alpha01];
        %i = 2;
            A(2,1:3) = [alpha10 1 alpha12];        
        %i = 3;
            A(3,2:4) = [alpha21 1 alpha23];        
        
        for i=4:NI-3
            A(i,i-1:i+1) = LHS;
        end
        
        %i = NI-2;
            A(NI-2,NI-3:NI-1) = [alpha23 1 alpha21];        
        %i = NI-1;
            A(NI-1,NI-2:NI)   = [alpha12 1 alpha10];        
        %i = NI;
            A(NI  ,NI-1:NI)   = [alpha01 1];        

A = inv(A);

for j=1:NJ
    for k=1:NK
        X = N(:,j,k);
        
        %i = 1;
            B(1)     = [a00 a01 a02 a03]*X(1:4)/dx^2;
        %i = 2;
            B(2)     = [a10 a11 a12 a13 a14]*X(1:5)/dx^2;
        %i = 3;
            B(3)     = [a20 a21 a22 a23 a24 a25]*X(1:6)/dx^2;
        
        for i=4:NI-3
            B(i)         = RHS*X(i-3:i+3);
        end
        
        %i = NI-2;
            B(NI-2)           = [a25 a24 a23 a22 a21 a20]*X(NI-5:NI)/(-dx)^2;
        %i = NI-1;
            B(NI-1)           = [a14 a13 a12 a11 a10]*X(NI-4:NI)/(-dx)^2;
        %i = NI;
            B(NI)             = [a03 a02 a01 a00]*X(NI-3:NI)/(-dx)^2;
            
        N(:,j,k) = A*B;
    end
end



N = ipermute(N, [dir, 1:dir-1, dir+1:ndim_N]);
