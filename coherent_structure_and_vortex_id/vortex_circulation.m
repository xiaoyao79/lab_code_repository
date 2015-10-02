function [Lcirc,Acirc]=vortex_circulation(U,V,x,y,pix,spacing,W_surface)
% The function vortex_circulation.m is written to solve for circulation of
% identified coherent structures based on the Line Integral and Area
% Integral approaches.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code:
%   1. Will compute circulation from velocity around the closed contour
%   2. Will compute circulation from the vorticity in the contour area
% Required inputs:
%   1. Velocity field in physcial units (U,V)
%   2. Vortex perimeter points (x,y; output from vortex_id.m)
%   3. pix (pixel magnification in meters/pixel)
%   4. spacing (spacing of PIV vectors (usually 4 or 8))
%
% Note: This code was modified from v1.0 by Chris Weiland 10.02.2007
% 
% v1: Brett Meyers, 11/21/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the size of x vector from vortex_id.m
[~,q2]  = size(x);
[s1,s2] = size(U);

% interpolate velocity fields
[xmat,ymat]=meshgrid(1:s2,1:s1);
Vx_int=interp2(xmat,ymat,U(:,:,1),x(1:q2-1),y(1:q2-1));
Vy_int=interp2(xmat,ymat,V(:,:,1),x(1:q2-1),y(1:q2-1));

% velocity vector information
Vmag=(Vx_int.^2+Vy_int.^2).^(1/2); % velocity vector magnitude
Vang=atan2(Vy_int,Vx_int); % direction of velocity vector
od=find(Vang<0); % ensure we work in absolute coordinate system.. atan2 treats 330 deg as -30 deg
Vang(od)=Vang(od)+2*pi;

% get distance between successive points and angle for line integral calculation
dx=diff(x).*-1; % incremental x-position
dy=diff(y).*-1; % incremental y-position
ds=(pix*spacing*(dx.^2+dy.^2).^0.5); % incremental distance converted to meters
ds_ang=atan2(dy,dx); % angle between sucessive points
od=find(ds_ang<0); % ensure we work in absolute coordinate system.. atan2 treats 330 deg as -30 deg
ds_ang(od)=ds_ang(od)+2*pi;

% angle between ds and V vector
ang_diff=abs(ds_ang-Vang);
od=find(ang_diff<0);
ang_diff(od)=ang_diff(od)+2*pi;

% take sum(V*ds*cos(theta))...this is the line integral over the closed contour
Lcirc=sum(Vmag.*ds.*cos(ang_diff));

% Calculated Circulation using Surface integral
[W_loc_x, W_loc_y]=find(~isnan(W_surface));
for i=1:size(W_loc_x);
    W_sum(i)=W_surface(W_loc_x(i), W_loc_y(i))*pix^2;
end
Acirc=sum(W_sum);