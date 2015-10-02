function [N]=socdiff2_bc(N,dx,dir,bcmask)
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
    
if NI<4
        DN(1     ,:,:) = (N(3   ,:,:) -2*N(2     ,:,:) + N(1     ,:,:)) / (dx^2);
        DN(2:NI-1,:,:) = (N(3:NI,:,:) -2*N(2:NI-1,:,:) + N(1:NI-2,:,:)) / (dx^2);
        DN(NI    ,:,:) = (N(  NI,:,:) -2*N(  NI-1,:,:) + N(  NI-2,:,:)) / (dx^2);
else
        DN(1     ,:,:) = ( -N(4   ,:,:) +4*N(3     ,:,:) -5*N(2     ,:,:) +2*N(1   ,:,:)) / (dx^2);
        DN(2:NI-1,:,:) = (  N(3:NI,:,:) -2*N(2:NI-1,:,:) +  N(1:NI-2,:,:)               ) / (dx^2);
        DN(NI    ,:,:) = (2*N(  NI,:,:) -5*N(  NI-1,:,:) +4*N(  NI-2,:,:) -  N(NI-3,:,:)) / (dx^2);    
end
    
%keyboard
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
    
    DN = zeros(NI,NJ,NK);
    
    %2o central diff
    [I,J] = find(bcmask & P1 & M1);
    DN(s2i([NI,NJ*NK],I,J)) = ( N(s2i([NI,NJ*NK],I+1,J)) - 2*N(s2i([NI,NJ*NK],I-1,J)) + N(s2i([NI,NJ*NK],I-1,J)) ) / (dx^2);
    %keyboard
    %2o forward diff (is it 2o?)
    [I,J] = find(bcmask & P1 & P2 & ~M1);
    DN(s2i([NI,NJ*NK],I,J)) = ( N(s2i([NI,NJ*NK],I+2,J)) - 2*N(s2i([NI,NJ*NK],I+1,J)) + N(s2i([NI,NJ*NK],I,J)) ) / (dx^2);
    %2o backward diff (is it 2o?)
    [I,J] = find(bcmask & M1 & M2 & ~P1);
    DN(s2i([NI,NJ*NK],I,J)) = ( N(s2i([NI,NJ*NK],I-2,J)) - 2*N(s2i([NI,NJ*NK],I-1,J)) + N(s2i([NI,NJ*NK],I,J)) ) / (dx^2);
    %1o forward diff -- could make assumptions, assume 0 to match 1st deriv
    [I,J] = find(bcmask & P1 & ~P2 & ~M1);
    %DN(s2i([NI,NJ*NK],I,J)) = 0;
    %1o backward diff - could make assumptions, assume 0 to match 1st deriv
    [I,J] = find(bcmask & M1 & ~M2 & ~P1);
    %DN(s2i([NI,NJ*NK],I,J)) = 0;
    %singletons and bcmask == 0 will be DN=0
end
N = ipermute(DN, [dir, 1:dir-1, dir+1:ndim_N]);

function N=s2i(siz,I,J)
%faster version of sub2ind
%converts I,J subscript into linear index into array
N=siz(1)*(J-1)+I;