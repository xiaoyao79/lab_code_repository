function [N]=socdiff_bc(N,dx,dir,bcmask)
% updated form of socdiff that only returns differences within bcmask, 0
% outside.

%sam found error 2008/05/14, k loop should be commented

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

if nargin < 4 %no bcmask defined
    
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
    
    
    if NI < 3
        error('stencil size is 3 pts in active direction')
    end
    
    DN = zeros(NI,NJ,NK);
    
    DN(1     ,:,:) = ( -N( 3   ,:,:) +4*N(  2 ,:,:) -3*N(1     ,:,:)) / (2*dx);
    DN(2:NI-1,:,:) = (  N( 3:NI,:,:)                -  N(1:NI-2,:,:)) / (2*dx);
    DN(NI    ,:,:) = (  N( NI-2,:,:) -4*N(NI-1,:,:) +3*N(NI    ,:,:)) / (2*dx);
    %DN(NI    ,:,:) = (3*N(NI   ,:,:) -4*N(NI-1,:,:) +  N(  NI-2,:,:)) / (2*dx);

else %use bcmask to limit area finite differenced
    %disp('using bcmask')
    if size(bcmask) ~= size(N)
        error('bcmask must be the same size as N')
    end
    
    N = permute(N, [dir, 1:dir-1, dir+1:ndim_N]);
    size_N = size(N);   %find again for permuted matrix
    ndim_N = ndims(N);
    
    bcmask = permute(bcmask, [dir, 1:dir-1, dir+1:ndim_N]);
    
    NI = size_N(1);
    NJ = size_N(2);
    if ndim_N == 3
        NK = size_N(3);
    else
        NK = 1;
    end
    
    if NI < 3
        error('stencil size is 3 pts in active direction')
    end
  
    %setup matrices to test existence of neighbors
    P1 = zeros(size_N);
    P2 = zeros(size_N);
    M1 = zeros(size_N);
    M2 = zeros(size_N);
    
    P1(1:end-1,:)     = bcmask(2:end,:); 
    P2(1:end-2,:)     = bcmask(3:end,:);
    M1(2:end,:)       = bcmask(1:end-1,:);
    M2(3:end,:)       = bcmask(1:end-2,:);
    %P1(end,:)         = 0;
    %P2(end-1:end,:)   = 0;
    %M1(1,:)           = 0;
    %M2(1:2,:)         = 0;
    
    DN = zeros(NI,NJ,NK);
    
%     keyboard
    
    %2o central diff
    [I,J] = find(bcmask & P1 & M1);
    DN(s2i([NI,NJ*NK],I,J)) = ( N(s2i([NI,NJ*NK],I+1,J)) - N(s2i([NI,NJ*NK],I-1,J)) ) / (2*dx);
    %keyboard
    %2o forward diff
    [I,J] = find(bcmask & P1 & P2 & ~M1);
    DN(s2i([NI,NJ*NK],I,J)) = ( -N(s2i([NI,NJ*NK],I+2,J)) + 4*N(s2i([NI,NJ*NK],I+1,J)) - 3*N(s2i([NI,NJ*NK],I,J)) ) / (2*dx);
    %2o backward diff
    [I,J] = find(bcmask & M1 & M2 & ~P1);
    DN(s2i([NI,NJ*NK],I,J)) = (  N(s2i([NI,NJ*NK],I-2,J)) - 4*N(s2i([NI,NJ*NK],I-1,J)) + 3*N(s2i([NI,NJ*NK],I,J)) ) / (2*dx);
    %1o forward diff
    [I,J] = find(bcmask & P1 & ~P2 & ~M1);
    DN(s2i([NI,NJ*NK],I,J)) = ( N(s2i([NI,NJ*NK],I+1,J)) - N(s2i([NI,NJ*NK],I,J)) ) / (dx);
    %1o backward diff
    [I,J] = find(bcmask & M1 & ~M2 & ~P1);
    DN(s2i([NI,NJ*NK],I,J)) = ( N(s2i([NI,NJ*NK],I,J)) - N(s2i([NI,NJ*NK],I-1,J)) ) / (dx);
    %singletons and bcmask == 0 will be DN=0
    
end

%keyboard

N = ipermute(DN, [dir, 1:dir-1, dir+1:ndim_N]);

function N=s2i(siz,I,J)
%faster version of sub2ind
%converts I,J subscript into linear index into array
N=siz(1)*(J-1)+I;

