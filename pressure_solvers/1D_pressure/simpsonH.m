function [Q] = simpsonH(F,h);
%[Q] = function SimpsonH(F,h);
%   Perform numerical integration of function f(x) on x=[a,b] sampled at 
%   k evenly spaced grid points on interval h=(b-a)/(k-1).
%
%   If k is odd, normal composite Simpson's rule will be used.
%   If k is even, Simpson's 3/8's rule will be used to complete the final
%   interval.
%
%   Integration will be performed on the columns of F.

[M,N] = size(F);

%could add check to transpose if data is in single row instead of column
if M <3
    error('F must have at least 3 points to integrate')
end

W = zeros(1,M);

if M==4             %use only 3/8's rule
    W = [1 3 3 1];
    Q = 3/8*h*W*F;
elseif mod(M,2)     %M is odd - use only Simpson's rule
    W(1) = 1;
    for m=2:2:M-1
        W(m) = 4;    
    end
    for m=3:2:M-2
        W(m) = 2;
    end
    W(M) = 1;
    Q = h/3*W*F;
else                %M is even - need 3/8's rule to finish interval
    W(1) = 8;
    for m=2:2:M-4
        W(m) = 32;
    end
    for m=3:2:M-5
        W(m) = 16;
    end
    W(M-3)=17;
    W(M-2)=27;
    W(M-1)=27;
    W(M) = 9;
    Q = h/24*W*F;
end

