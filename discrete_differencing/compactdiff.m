function [N]=compactdiff(N,dx,dir)

% if nargin<3
%     %from Lele's original work
%     alpha = 1/3; beta=0;    a=14/9;  b=1/9;     c=0;       %tridiagonal,5pt stencil,6th order
%     alpha = 3/8; beta=0;    a=25/16; b=1/5;     c=-1/80;   %tridiagonal,7pt stencil,8th order
%     alpha = 1/2; beta=1/20; a=17/12; b=101/150; c=1/100;   %10th order
% end

%Kim&Lee 1996 wave-optimized compact schemes
%triadiagonal, 7pt stencil
    %a=1.545790417; b=0.434249728; c=-0.078236437; alpha=0.450901855; beta=0;    %2nd order
    %a=1.551941906; b=0.361328195; c=-0.042907397; alpha=0.435181352; beta=0;    %4th order
    a=1.568098212; b=0.271657107; c=-0.022576781; alpha=0.408589269; beta=0;    %6th order
%pentadiagonal, 7pt stencil
    %a=1.265667929; b=1.079904285; c=0.053798648; alpha=0.596631925; beta=0.103053504; %2nd order
    %a=1.280440844; b=1.049309076; c=0.044465832; alpha=0.589595521; beta=0.097512355; %4th order
    %a=1.323482375; b=0.944394243; c=0.027596356; alpha=0.566458285; beta=0.081278202; %6th order
    %a=1.373189728; b=0.814447053; c=0.016707870; alpha=0.537265947; beta=0.064906379; %8th order
%tridiagonal 6th order seems to be best performing by their paper

%BC - Tridiagonal
%i=0, 2nd order
a00 = -2.67344438910815;  a01 = 1.46806676496733;
a02 =  1.38268870248505;  a03 = -0.177311078344225;
alpha01 = 2.70151093490474;
%i=1, 4th order
a10 = -0.508867575457385; a11 = -0.702987853336675;
a12 =  1.04038536544838; a13 = 0.186747203650676;
a14 = -0.015277140304991;  
alpha10 = 0.153204878183875; alpha12 = 0.723711049108264;
%i=2, 6th order
a20 = -0.013127263621621; a21 = -0.603802922173413;
a22 = -0.439515424684709; a23 =  0.96090920472974;
a24 =  0.101030348558563; a25 = -0.0054939428085588;
alpha21 = 0.223454477162156; alpha23 = 0.553091045675688;


%%leftover coefficients for pentadiagonal schemes - organize based on paper
% 2. Pentadiagonal Schemes:
% i = 0: Second Order:
% ow = - 2.95516745786296, 00,1'
% 0^ = 4.28093227034817, ow
% :
% 00,1=4.57321732196853, fa.
% i=l: Fourth Order:
% aijo = - 0.643755519081585, 01,1 =
% 01,2 = 1.39308006947385, op:
% 01,4 = -0.055981044934069, 01,0
% 01,2 = 0.046406522760991, fo
% 
% 
% -1.63175038219495
% 0.305985569709741
% 2.27485354566209
% -0.215562412498565
% -0.477781092959631
% -0.204356208611126
% -0.337432463538152

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

RHS = [ -c/6 -b/4 -a/2 0 a/2 b/4 c/6 ]/dx;
%LHS = [ beta alpha 1 alpha beta ]; %pentadiagonal
LHS = [ alpha 1 alpha ];  %tridiagonal

%setup A matrix:
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
%keyboard
           
%pre-solve inverse
%A = A\speye(NI,NI);
A = inv(A);
% Ad = spdiags(A);
% BB = Ad(:,3);
% DD = Ad(:,2);
% AA = Ad(:,1);

for j=1:NJ
    for k=1:NK
        X = N(:,j,k);
        
        %i = 1;
            B(1)     = [a00 a01 a02 a03]*X(1:4)/dx;
        %i = 2;
            B(2)     = [a10 a11 a12 a13 a14]*X(1:5)/dx;
        %i = 3;
            B(3)     = [a20 a21 a22 a23 a24 a25]*X(1:6)/dx;
        
        for i=4:NI-3
            B(i)         = RHS*X(i-3:i+3);
        end
        
        %i = NI-2;
            B(NI-2)           = [a25 a24 a23 a22 a21 a20]*X(NI-5:NI)/(-dx);
        %i = NI-1;
            B(NI-1)           = [a14 a13 a12 a11 a10]*X(NI-4:NI)/(-dx);
        %i = NI;
            B(NI)             = [a03 a02 a01 a00]*X(NI-3:NI)/(-dx);
            
        N(:,j,k) = A*B;                                    %much slower
        %N(:,j,k) = TDS1(1,NI,Ad(:,3),Ad(:,2),Ad(:,1),B);   %almost as fast
        %N(:,j,k) = TDS1(1,NI,BB,DD,AA,B);                   %seems to be fastest
    end
end

%keyboard

N = ipermute(N, [dir, 1:dir-1, dir+1:ndim_N]);
