function [V1,V2,D,A] = KLDs2v(U1,U2,GRIDTYPE,w_yj,c_j)
%[phi, lambda, a] = KLD2v(U1,U2,GRIDTYPE,W,C)
% Perform Karhunen-Loeve / Proper-Orthogonal decomposition on the two 
% dimensional vector signal U_i integrated on an associated GRIDTYPE.
% The gridtype can be 'even', 'cosine', or 'custom'.  For a custom grid, a
% weight function W and a quadrature function C must be entered, defined at
% the grid points of U.  W and C should be the same length as Ui.

% U to be decomposed needs to be reordered as a K by NT array, where 
% If U is originally 2D, reorder U as a single strip first.  
% 

if nargin<3
    GRIDTYPE='even';
end

[K1,K2,NT] = size(U1);
if size(U2)~=[K1,K2,NT]
    error('U1 and U2 must be the same sizes')
end

switch lower(GRIDTYPE(1:2))
    case {'ev','re','sq'}   %even, rectangular, square
        w_yj = ones(K1,K2);
        c1_j     = ones(K1,1);
        c1_j(1)  = 0.5;
        c1_j(K1) = 0.5;
        c2_j     = ones(K2,1);
        c2_j(1)  = 0.5;
        c2_j(K2) = 0.5;
        c_j = c1_j*c2_j.';
        clear c1_j c2_j;
    case {'ch','co'}        %chebyshev, cosine
        y1j = [cos(pi*[0:K1]/K1)].';
        w1_yj = sqrt(1-y1j.^2).^-1;
        y2j = [cos(pi*[0:K2]/K2)].';
        w2_yj = sqrt(1-y2j.^2).^-1;
        w_yj = w1_yj*w2_yj.';   

        c1_j     = pi/K1*ones(K1,1);
        c1_j(1)  = pi/2/K1;
        c1_j(K1) = pi/2/K1;
        c2_j     = pi/K2*ones(K2,1);
        c2_j(1)  = pi/2/K2;
        c2_j(K2) = pi/2/K2;
        c_j = c1_j*c2_j.';

        clear w1_yj w2_yj c1_j c2_j;
    case {'cu'}             %custom
        if nargin <5
            error('A weight function and quadrature coefficients must be defined for a ''custom'' GRIDTYPE')
        end

        if length(w_yj)~= [K1,K2] or length(c_j) ~= [K1,K2]
            error('The size of W and C must be the same as U')
        end
    otherwise
        error('GRIDTYPE must be ''even'', ''cosine'', or ''custom''')
end

%unfold 2D fields into 1D strips
U1    = reshape(permute(U1,[2,1,3]),K1*K2,NT,1);
U2    = reshape(permute(U2,[2,1,3]),K1*K2,NT,1);
w_yj = reshape(permute(w_yj,[2,1]),K1*K2,1);
c_j  = reshape(permute( c_j,[2,1]),K1*K2,1);

%concatenate fields
[V,D,A] = KLDs([U1;U2],'custom',[w_yj;w_yj],[c_j;c_j]);

%need to reorder 1-D eigenvectors into 2D eigenmodes
%U = permute(reshape(U,K2,K1,NT),[2,1,3]);
V1 = V(      1:  K1*K2,:);
V2 = V(K1*K2+1:2*K1*K2,:);
V1 = permute(reshape(V1,K2,K1,NT),[2,1,3]);
V2 = permute(reshape(V2,K2,K1,NT),[2,1,3]);
