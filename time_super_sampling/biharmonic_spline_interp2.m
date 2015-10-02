function ZI = biharmonic_spline_interp2(X,Y,Z,XI,YI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2D Biharmonic spline interpolation implemented from:
%
%   Sandwell, D. T. (1987), Biharmonic spline interpolation of GEOS-3 and
%       SEASAT altimeter data, Geophysical Research Letters, Vol. 2, 
%       p. 139 ? 142.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run an example if no input arguments are found
if nargin ~= 5
    fprintf('Running Peaks Example \n');
    X = rand(100,1)*6 - 3;
    Y = rand(100,1)*6 - 3;
    Z = peaks(X,Y);
    [XI,YI] = meshgrid(-3:.25:3);
end
%   Check to make sure sizes of input arguments are correct/consistent
if size(Z,1) < size(Z,2)
    error('X, Y and Z should all be column vectors !');
end
if (size(Z,1) ~= size(X,1)) || (size(Z,1) ~= size(Y,1))
    error('Length of X, Y, and Z must be equal !');
end
if (size(XI,1) ~= size(YI,1)) || (size(XI,2) ~= size(YI,2))
    error('Size of XI and YI must be equal !');
end
%   Initialize output
ZI = zeros(size(XI));
%   Compute GG matrix for GG*m = d inversion problem
GG = zeros(length(Z),length(Z));
for i = 1 : length(Z)
    for j = 1 : length(Z)
        if i ~= j
            magx = sqrt((X(i)-X(j))^2 + (Y(i)-Y(j))^2);
            if magx >= 1e-7
                GG(i,j) = (magx^2) * (log(magx)-1);
            end
        end
    end
end
%   Compute model "m" where data "d" is equal to "Z"
m = GG\Z;
%   Find 2D interpolated surface through irregular/regular X, Y grid points 
gg = zeros(size(m));
for i = 1 : size(ZI,1)
    for j = 1 : size(ZI,2)
        for k = 1 : length(Z)
            magx = sqrt((XI(i,j)-X(k))^2 + (YI(i,j)-Y(k))^2);
            if magx >= 1e-7
                gg(k) = (magx^2) * (log(magx)-1);
            end 
        end
        ZI(i,j) = sum(gg.*m);
    end
end
%   Plot result if running example or if no output arguments are found
if nargin ~= 5 || nargout ~= 1
    figure;
    subplot(3,1,1);
        [XE,YE] = peaks(-3:.25:3);
        ZE = peaks(XE,YE);
        mesh(XE,YE,ZE); hold on;
        scatter3(X,Y,Z,'filled'); hold off;
        title('Peaks');
        axis([-3 3 -3 3 -5 5]);
        caxis([-5 5]); colorbar;
    subplot(3,1,2);
        mesh(XI,YI,ZI); hold on;
        scatter3(X,Y,Z,'filled'); hold off;
        title('Peaks Interpolated'); 
        axis([-3 3 -3 3 -5 5]);   
        caxis([-5 5]); colorbar;
    subplot(3,1,3);
        mesh(XI,YI,ZI-ZE); hold on;
        scatter3(X,Y,Z-Z,'filled'); hold off;
        title('Peaks Interpolated Difference'); 
        axis([-3 3 -3 3 -5 5]);  
        caxis([-5 5]); colorbar; 
end