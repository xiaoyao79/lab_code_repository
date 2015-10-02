function [dudx ]=gradient_dudxdvdy_compact_rich(N,dx,dir)
% This function calculates gradients using the 4th order compact-richardson
% scheme introduced by A. Etebari and P. Vlachos in "Improvements on
% the accuracy of derivative estimation from DPIV velocity measurements"
% Experiments in Fluids (2005)

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
ADN = zeros(NI,NJ,NK);

A=[1239 272 1036 -69 0];
k=[1 2 4 8];

%% Parameters for Boundary formulation for the First
%% Derivative (Lele et al 1992)
alpha=1;
d1=0;
a1=-(3+alpha+2*d1)/2;
b1=2+3*d1;
c1=-(1-alpha+6*d1)/2;

for j=1:NJ
    current_count=1;
    for m=1:4;
        
        c1=0;
        while c1<=k
%             keyboard
            N_current(:,1)=N(k(m)+c1:k(m):size(N,1),j);
            
            a(:,1)=1/4*ones(1,size(N_current,1));
            a(1)=0;
            a(size(N_current,1))=0;
            b(:,1)=ones(1,size(N_current,1));
            b(1)=1;
            b(size(N_current,1))=1;
            c(:,1)=1/4*ones(1,size(N_current,1));
            c(1)=0;
            c(size(N_current,1))=0;
            for l=2:size(N_current,1)-1
                d(l,1)=1.5*(N_current(l+1,1)-N_current(l-1,1))/(2*k(m)*dx);
            end
            
            d(1,1)=1/(k(m)*dx)*(N_current(2,1)-N_current(1,1));
            d(size(N_current,1),1)=1/(k(m)*dx)*(N_current(size(N_current,1),1)-N_current(size(N_current,1)-1,1));
            
            [y]=tridiagSolve(a,b,c, d);   %% Solve tridiagonal Matrix
            ADN(k(m)+c1:k(m):size(N,1),j)=A(m+1)*y'+ADN(k(m)+c1:k(m):size(N,1),j); %% Calculate current sum of (A*U')
            clear a b c d y
            clear N_current
            c1=c1+1;
        end
        
    end
    DN=1/A(1)*ADN;
end
dudx = ipermute(DN, [dir, 1:dir-1, dir+1:ndim_N]);