function [V,D,A] = KLD(U,GRIDTYPE,w_yj,c_j)
%[phi, lambda, a] = KLD(U,GRIDTYPE,W,C)
% Perform Karhunen-Loeve / Proper-Orthogonal decomposition on the one 
% dimensional signal U, integrated on an associated GRIDTYPE.
% The gridtype can be 'even', 'cosine', or 'custom'.  For a custom grid, a
% weight function W and a quadrature function C must be entered, defined at
% the grid points of U.  W and C should be the same length as U.

% U to be decomposed needs to be reordered as a K by NT array, where 
% If U is originally 2D, reorder U as a single strip first.  
% 

if nargin<2
    GRIDTYPE='even';
end

[K,NT] = size(U);

switch lower(GRIDTYPE(1:2))
    case {'ev','re','sq'}   %even, rectangular, square
        w_yj = ones(K,1);
        c_j  = ones(K,1);
        c_j(1) = 0.5;
        c_j(K) = 0.5;
    case {'ch','co'}        %chebyshev, cosine
        yj = [cos(pi*[0:K]/K)].';
        w_yj = sqrt(1-yj.^2).^-1;
        c_j  = pi/K*ones(K,1);
        c_j(1) = pi/2/K;
        c_j(K) = pi/2/K;
    case {'cu'}             %custom
        if nargin <4
            error('A weight function and quadrature coefficients must be defined for a ''custom'' GRIDTYPE')
        end
        w_yj = w_yj(:); %make vertical
        c_j  = c_j(:);

        if length(w_yj)~= K or length(c_j) ~=K
            error('The length of W and C must be the same as U')
        end
    otherwise
        error('GRIDTYPE must be ''even'', ''cosine'', or ''custom''')
end

UUt_avg = zeros(K,K);

for n=1:NT
    UUt_avg = UUt_avg + U(:,n)*U(:,n)';
end

UUt_avg = UUt_avg/NT;

sw = sqrt(w_yj);
sc = sqrt(c_j);

% %scaling factors needed to keep quadrature matrix Hermetian
% ww = sw*sw.';
% cc = sc*sc.';
% KK = UUt_avg.*cc./ww; 
% [V,D] = eig(KK);

% KK = UUt_avg.*(sc*sc.')./(sw*sw.'); 
% [V,D] = eig(KK);

[V,D] = eig(UUt_avg.*(sc*sc.')./(sw*sw.'));

%a little worried about whether this holds up generally
if isreal(UUt_avg)
    V = real(V);
    D = real(D);
end

%reverse scaling to natural magnitudes
for n=1:K
    V(:,n) = V(:,n).*sw./sc;
end

%keyboard
[D,Di] = sort(diag(D),1,'descend');
V = V(:,Di);

% % had to remove this - either needed to store original V, or add some more
% % if statements
% %only calculate time coefficients if we were asked for them
% if nargout<3
%     return
% end

%a_t,k = sum_j ( U_t,j * V*_j,k/w_j ) = U_t,j*V*_j,k/w_j
%a_k,t = (V*_j,k/w_j)T * U_j,t = V*_k,j/w_j * U_j,t
A = (V.*repmat(c_j./w_yj,[1,K]))' * U;

%correct for normalization on V (include weights??)
Vnorm = sqrt(sum((V).^2,1));
for i=1:NT
    V(:,i) = V(:,i)/Vnorm(i);
    A(i,:) = A(i,:)*Vnorm(i);
    %D(i) = D(i)*Vnorm(i);
end




