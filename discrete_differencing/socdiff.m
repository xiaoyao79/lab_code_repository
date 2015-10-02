function [N]=socdiff(N,dx,dir)
%

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
%for j=1:NJ
%    for k=1:NK
        DN(1     ,:,:) = ( -N( 3   ,:,:) +4*N(  2 ,:,:) -3*N(1     ,:,:)) / (2*dx);
        DN(2:NI-1,:,:) = (  N( 3:NI,:,:)                -  N(1:NI-2,:,:)) / (2*dx);
        DN(NI    ,:,:) = (3*N(NI   ,:,:) -4*N(NI-1,:,:) +  N(  NI-2,:,:)) / (2*dx);
%    end
%end

%keyboard

N = ipermute(DN, [dir, 1:dir-1, dir+1:ndim_N]);
