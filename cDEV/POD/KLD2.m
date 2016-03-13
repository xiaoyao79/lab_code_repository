function [V,D,A] = KLD2(U,GRIDTYPE,w_yj,c_j)
%[phi, lambda, a] = KLD2(U,GRIDTYPE,W,C)
% Perform Karhunen-Loeve / Proper-Orthogonal decomposition on the two 
% dimensional signal U, integrated on an associated GRIDTYPE.
% The gridtype can be 'even', 'cosine', or 'custom'.  For a custom grid, a
% weight function W and a quadrature function C must be entered, defined at
% the grid points of U.  W and C should be the same length as U.

% U to be decomposed needs to be reordered as a K by NT array, where 
% If U is originally 2D, reorder U as a single strip first.  


% load Case_5_C_u_vect.mat;

if nargin<2
    GRIDTYPE='even';
end

[K1,K2,NT] = size(U);

switch lower(GRIDTYPE(1:2))
    case {'ev','re','sq'}   %even, rectangular, square
        w_yj = ones(K1,K2);
        c1_j     = ones(K1);
        c1_j(1)  = 0.5;
        c1_j(K1) = 0.5;
        c2_j     = ones(K2);
        c2_j(1)  = 0.5;
        c2_j(K2) = 0.5;
        c_j = c1_j*c2_j';
        clear c1_j c2_j;
    case {'ch','co'}        %chebyshev, cosine
        y1j = [cos(pi*[0:K1]/K1)].';
        w1_yj = sqrt(1-y1j.^2).^-1;
        y2j = [cos(pi*[0:K2]/K2)].';
        w2_yj = sqrt(1-y2j.^2).^-1;
        w_yj = w1_yj*w2_yj.';   

        c1_j     = pi/K*ones(K1);
        c1_j(1)  = pi/2/K;
        c1_j(K1) = pi/2/K;
        c2_j     = pi/K*ones(K2);
        c2_j(1)  = pi/2/K;
        c2_j(K2) = pi/2/K;
        c_j = c1_j*c2_j.';

        clear w1_yj w2_yj c1_j c2_j;
    case {'cu'}             %custom
        if nargin <4
            error('A weight function and quadrature coefficients must be defined for a ''custom'' GRIDTYPE')
        end

        if length(w_yj)~= [K1,K2] or length(c_j) ~= [K1,K2]
            error('The size of W and C must be the same as U')
        end
    otherwise
        error('GRIDTYPE must be ''even'', ''cosine'', or ''custom''')
end

%unfold 2D fields into 1D strips
U    = reshape(permute(U,[2,1,3]),K1*K2,NT,1);
w_yj = reshape(permute(w_yj,[2,1]),K1*K2,1);
c_j  = reshape(permute( c_j,[2,1]),K1*K2,1);

[V,D,A] = KLD(U,'cosine',w_yj,c_j);

%need to reorder 1-D eigenvectors into 2D eigenmodes
%U = permute(reshape(U,K2,K1,NT),[2,1,3]);
V = permute(reshape(V,K2,K1,K1*K2),[2,1,3]);
