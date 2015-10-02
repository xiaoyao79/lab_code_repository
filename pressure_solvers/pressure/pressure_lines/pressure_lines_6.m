function [p]=pressure_lines(x,y,u,v,Re,dt,solver,bclist,Pref,Px,Py)
%function [p]=pressure_solve(x,y,u,v,Re,dt,[solver],[bclist],[Pref, Px, Py])
%
% Returns p(x,y,t) from u(x,y,t) and v(x,y,t) using a multi-directional 
% Guass-Seidel style of line integration on pressure terms of the Navier-\
% Stokes Equations.  By default, a conservative formulation.
%
% SOLVER must be either "conservative" or "standard". The default is a
% conservative formulation for the convective terms
%   "standard" Form:
%     dP/dx = 1/Re * (d2u/dx2+d2u/dy2) - ( du/dt + u*du/dx + v*du/dy )
%     dP/dy = 1/Re * (d2v/dx2+d2v/dy2) - ( dv/dt + u*dv/dx + v*dv/dy )
%
%   "conservative" Form:
%     dP/dx = 1/Re * (d2v/dx2+d2v/dy2) - ( du/dt + d/dx[uu] + d/dy[uv] + u*dw/dz )
%     dP/dy = 1/Re * (d2v/dx2+d2v/dy2) - ( dv/dt + d/dx[vu] + d/dy[vv] + v*dw/dz )
%     dw/dz = -( du/dx + dv/dy )
%
% The solver assumes that X and Y are uniformly spaced, but dX and dY do 
% not have to be the same, and are determined from X and Y
%
% All variables should be normalized by:
%   u = U'/U0, v = V'/U0, x = X'/L, y = Y'/L
%   t = T'/T , T=L/U0, Re = L*U0/nu
%   p = P'/P0, P0 = rho*U0^2
%
% BCLIST must have the following form:
% { [WALL1],[WALL2], ... }
% where each of the [WALLS] is a N*2 list of Neumann BC that has the
% form [X1,Y1; X2,Y2; ... ; Xn,Yn].  Each X and Y should have a value 
% corresponding to an indice from the matrices U and V.  The lengths 
% N need not be equal for each wall.  The list of points should be adjacent
% and closed (though the last point does not have match the first).  The
% fluid boundary is the right hand side of each wall, i.e. the flow lies
% inside a clockwise loop, and outside a counterclockwise one.
%
% If BCLIST is not defined, the default is to use the domain walls as BC's
%
% PREF is a normalized pressure measured at index pt [PX,PY]
% If PREF is a single value, it is assumed that the pressure at [PX,PY] 
% is constant for all time steps.  Otherwise it must have the same number 
% of time steps as U and V. PREF is optional, but if not defined, the
% pressure will float between time steps and cannot be compared.

% converted derivatives to compact schemes - this breaks BC checks, so
%  removed them
% updated derivative calc on boundaries of domains...
% updated derivative calc on given boundary conditions
% added pressure Dirichlet BC Pref in interior at Px and Py
% added switch for integration equations, removed ideal vs viscous
% might consider making viscous terms conservative too
% carried correction on d2vdxdy from pressure_direct2.m
% added check to see if enough time steps (>=7) for compactdiff, else go
%  back to central differences

% try

%figure(1234),quiver(x,y,u(:,:,1),v(:,:,1)),pause

NI = size(x,1);
NJ = size(y,2);
NT = size(u,3);

%dynamic boundary conditions based on function arguments
%Neumann-type pressure gradient BC's
if nargin>=7
    if nargin>=8
        if nargin>=9
            if nargin<11
                error('If PREF is defined, PX and PY must also be given');
            end
            %disp('PREF, PX, and PY are defined')
            Pref = Pref(:); %trick to make Pref vertical
            if (length(Pref)==1)
                Pref = Pref*ones(size(u,3),1);
            else
                if length(Pref)~=size(u,3)
                    error('PREF must have length 1 or equal to number of flow fields in U and V');
                end
            end
        end
        %disp('BCLIST is defined')
        if iscell(bclist)
            %disp('BCLIST is a cell')
            if size(bclist) == 0;
                disp('BCLIST must have at least 1 element')  
                error('BCLIST must have the following form: { [x1,y1;...],[x1,y2;...], ... }')
            end 
        else %can't use this BCLIST
            disp('BCLIST is not a cell')                        
            error('BCLIST must have the following form: { [x1,y1;...],[x1,y2;...], ... }')
        end
    else 
        disp('BCLIST is not defined, using default boundary conditions')
            %define BCLIST to use the domain edges as BC's
            bclist{1} = [[      (1:NI-1)' , NJ*ones(NI-1,1) ];...    %north
                         [NI*ones(NJ-1,1) ,       (NJ:-1:2)'];...    %east  
                         [      (NI:-1:2)',    ones(NI-1,1) ];...    %south
                         [   ones(NJ-1,1),         (1:NJ-1)']];     %west
    end
else 
    %flowtype = 'viscous'
    solver = 'conservative'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
tic
restarts        = 20;
max_iterations  = 50;
tol             = 1e-5;

dx = x(2,1) - x(1,1);
dy = y(1,2) - y(1,1);

% x=0:dx:1;
% y=0:dy:1;
% NI=length(x)
% NJ=length(y)

%moved to top
% NI = size(x,1);
% NJ = size(y,2);
% NT = size(u,3);

% storage for source term
S = zeros(NI,NJ,NT);
% phi(i,j,k) where i is x-index, j is y-index, and k=1,2,3 for time-stepping
%discretization defined on the following neighboring points
AA = sparse(NI*NJ,NI*NJ);
An = zeros(NI,NJ);  
As = zeros(NI,NJ);
Aw = zeros(NI,NJ);
Ae = zeros(NI,NJ);
Ap = zeros(NI,NJ);
%store pressure in phi during solution
phi = zeros(NI*NJ,1);
%source term plus LHS terms not used in tridiagonal solution
% RHSi = zeros(NI,1);
% RHSj = zeros(NJ,1);
%storage for final pressure solutions
p = zeros(NI,NJ,NT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the source term 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_size = size(u);   %u and v should be the same size
M = u_size(1);      %number of points in x
N = u_size(2);      %number of points in y
T = u_size(3);      %number of points in t

%p = p(x,y,t)
%dp/dx = rho*[du/dt + u*du/dx + v*du/dy] = rho*f(x,y,t)
%dp/dy = rho*[dv/dt + u*dv/dx + v*dv/dy] = rho*g(x,y,t)
%dp/dx +dp/dy = rho*[f+g] = H(x,y,t)

%this is now defined above, derived from input data
% dx = x(2,1)-x(1,1);
% dy = y(1,2)-y(1,1);
% dt = 1;                 %should do something to enter this

%use 2nd order C-D schemes in all variables
%could replace these with implicit derivatives
%also, does not respect BC, should  fix that. (second pass below)
% dudt(:,:,1    ) = (u(:,:,2  ) - u(:,:,1    )) /    dt ;
% dudt(:,:,2:T-1) = (u(:,:,3:T) - u(:,:,1:T-2)) / (2*dt);
% dudt(:,:,    T) = (u(:,:,  T) - u(:,:,  T-1)) /    dt ;
% 
% dvdt(:,:,1    ) = (v(:,:,2  ) - v(:,:,1    )) /    dt ;
% dvdt(:,:,2:T-1) = (v(:,:,3:T) - v(:,:,1:T-2)) / (2*dt);
% dvdt(:,:,    T) = (v(:,:,  T) - v(:,:,  T-1)) /    dt ;

uu = u.*u;
uv = u.*v;
vv = v.*v;

%keyboard
if NT<7
    dudt(:,:,1    ) = ( -u(:,:,3  ) +4*u(:,:, 2 ) -3*u(:,:,1    )) / (2*dt);
    dudt(:,:,2:T-1) = (  u(:,:,3:T)               -  u(:,:,1:T-2)) / (2*dt);
    dudt(:,:,    T) = (3*u(:,:,  T) -4*u(:,:,T-1) +  u(:,:,  T-2)) / (2*dt);
    
    dvdt(:,:,1    ) = ( -v(:,:,3  ) +4*v(:,:, 2 ) -3*v(:,:,1    )) / (2*dt);
    dvdt(:,:,2:T-1) = (  v(:,:,3:T)               -  v(:,:,1:T-2)) / (2*dt);
    dvdt(:,:,    T) = (3*v(:,:,  T) -4*v(:,:,T-1) +  v(:,:,  T-2)) / (2*dt);
else
    dudt = compactdiff(u,dt,3);
    dvdt = compactdiff(v,dt,3);
end

dudx = compactdiff(u,dx,1);
% dudx(1    ,:,:) = ( -u(3  ,:,:) +4*u( 2 ,:,:) -3*u(1    ,:,:)) / (2*dx);
% dudx(2:M-1,:,:) = (  u(3:M,:,:)               -  u(1:M-2,:,:)) / (2*dx);
% dudx(M    ,:,:) = (3*u(M  ,:,:) -4*u(M-1,:,:) +  u(  M-2,:,:)) / (2*dx);

dvdx = compactdiff(v,dx,1);
% dvdx(1    ,:,:) = ( -v(3  ,:,:) +4*v( 2 ,:,:) -3*v(1    ,:,:)) / (2*dx);
% dvdx(2:M-1,:,:) = (  v(3:M,:,:)               -  v(1:M-2,:,:)) / (2*dx);
% dvdx(M    ,:,:) = (3*v(M  ,:,:) -4*v(M-1,:,:) +  v(  M-2,:,:)) / (2*dx);

dudy = compactdiff(u,dy,2);
% dudy(:,1    ,:) = ( -u(:,3  ,:) +4*u(:, 2 ,:) -3*u(:,1    ,:)) / (2*dy);
% dudy(:,2:N-1,:) = (  u(:,3:N,:)               -  u(:,1:N-2,:)) / (2*dy);
% dudy(:,    N,:) = (3*u(:,  N,:) -4*u(:,N-1,:) +  u(:,  N-2,:)) / (2*dy);

dvdy = compactdiff(v,dy,2);
% dvdy(:,1    ,:) = ( -v(:,3  ,:) +4*v(:, 2 ,:) -3*v(:,1    ,:)) / (2*dy);
% dvdy(:,2:N-1,:) = (  v(:,3:N,:)               -  v(:,1:N-2,:)) / (2*dy);
% dvdy(:,    N,:) = (3*v(:,  N,:) -4*v(:,N-1,:) +  v(:,  N-2,:)) / (2*dy);

d2udx2 = compactdiff2(u,dx,1);
% d2udx2(1    ,:,:) = (u(3  ,:,:) -2*u(2    ,:,:) + u(1    ,:,:)) / (dx^2);
% d2udx2(2:M-1,:,:) = (u(3:M,:,:) -2*u(2:M-1,:,:) + u(1:M-2,:,:)) / (dx^2);
% d2udx2(M    ,:,:) = (u(  M,:,:) -2*u(  M-1,:,:) + u(  M-2,:,:)) / (dx^2);

d2udy2 = compactdiff2(u,dy,2);
% d2udy2(:,1    ,:) = (u(:,3  ,:) -2*u(:,2    ,:) + u(:,1    ,:)) / (dy^2);
% d2udy2(:,2:N-1,:) = (u(:,3:N,:) -2*u(:,2:N-1,:) + u(:,1:N-2,:)) / (dy^2);
% d2udy2(:,    N,:) = (u(:,  N,:) -2*u(:,  N-1,:) + u(:,  N-2,:)) / (dy^2);

d2vdx2 = compactdiff2(v,dx,1);
% d2vdx2(1    ,:,:) = (v(3  ,:,:) -2*v(2    ,:,:) + v(1    ,:,:)) / (dx^2);
% d2vdx2(2:M-1,:,:) = (v(3:M,:,:) -2*v(2:M-1,:,:) + v(1:M-2,:,:)) / (dx^2);
% d2vdx2(M    ,:,:) = (v(M  ,:,:) -2*v(  M-1,:,:) + v(  M-2,:,:)) / (dx^2);

d2vdy2 = compactdiff2(v,dy,2);
% d2vdy2(:,1    ,:) = (v(:,3  ,:) -2*v(:,2    ,:) + v(:,1    ,:)) / (dy^2);
% d2vdy2(:,2:N-1,:) = (v(:,3:N,:) -2*v(:,2:N-1,:) + v(:,1:N-2,:)) / (dy^2);
% d2vdy2(:,    N,:) = (v(:,  N,:) -2*v(:,  N-1,:) + v(:,  N-2,:)) / (dy^2);

% d3udx2dx = compactdiff(d2udx2,dx,1);
% % d3udx2dx(1    ,:,:) = ( -d2udx2(3  ,:,:) +4*d2udx2( 2 ,:,:) -3*d2udx2(1    ,:,:)) / (2*dx);
% % d3udx2dx(2:M-1,:,:) = (  d2udx2(3:M,:,:)                    -  d2udx2(1:M-2,:,:)) / (2*dx);
% % d3udx2dx(M    ,:,:) = (3*d2udx2(M  ,:,:) -4*d2udx2(M-1,:,:) +  d2udx2(  M-2,:,:)) / (2*dx);
% 
% d3udy2dx = compactdiff(d2udy2,dx,1);
% % d3udy2dx(1    ,:,:) = ( -d2udy2(3  ,:,:) +4*d2udy2( 2 ,:,:) -3*d2udy2(1    ,:,:)) / (2*dx);
% % d3udy2dx(2:M-1,:,:) = (  d2udy2(3:M,:,:)                    -  d2udy2(1:M-2,:,:)) / (2*dx);
% % d3udy2dx(M    ,:,:) = (3*d2udy2(M  ,:,:) -4*d2udy2(M-1,:,:) +  d2udy2(  M-2,:,:)) / (2*dx);
% 
% d3vdx2dy = compactdiff(d2vdx2,dy,2);
% % d3vdx2dy(:,1    ,:) = ( -d2vdx2(:,3  ,:) +4*d2vdx2(:, 2 ,:) -3*d2vdx2(:,1    ,:)) / (2*dy);
% % d3vdx2dy(:,2:N-1,:) = (  d2vdx2(:,3:N,:)                    -  d2vdx2(:,1:N-2,:)) / (2*dy);
% % d3vdx2dy(:,    N,:) = (3*d2vdx2(:,  N,:) -4*d2vdx2(:,N-1,:) +  d2vdx2(:,  N-2,:)) / (2*dy);
% 
% d3vdy2dy = compactdiff(d2vdy2,dy,2);
% % d3vdy2dy(:,1    ,:) = ( -d2vdy2(:,3  ,:) +4*d2vdy2(:, 2 ,:) -3*d2vdy2(:,1    ,:)) / (2*dy);
% % d3vdy2dy(:,2:N-1,:) = (  d2vdy2(:,3:N,:)                    -  d2vdy2(:,1:N-2,:)) / (2*dy);
% % d3vdy2dy(:,    N,:) = (3*d2vdy2(:,  N,:) -4*d2vdy2(:,N-1,:) +  d2vdy2(:,  N-2,:)) / (2*dy);

duudx = compactdiff(uu,dx,1);
% duudx(1    ,:,:) = ( -uu(3  ,:,:) +4*uu( 2 ,:,:) -3*uu(1    ,:,:)) / (2*dx);
% duudx(2:M-1,:,:) = (  uu(3:M,:,:)                -  uu(1:M-2,:,:)) / (2*dx);
% duudx(M    ,:,:) = (3*uu(M  ,:,:) -4*uu(M-1,:,:) +  uu(  M-2,:,:)) / (2*dx);

duvdx = compactdiff(uv,dx,1);
% duvdx(1    ,:,:) = ( -uv(3  ,:,:) +4*uv( 2 ,:,:) -3*uv(1    ,:,:)) / (2*dx);
% duvdx(2:M-1,:,:) = (  uv(3:M,:,:)                -  uv(1:M-2,:,:)) / (2*dx);
% duvdx(M    ,:,:) = (3*uv(M  ,:,:) -4*uv(M-1,:,:) +  uv(  M-2,:,:)) / (2*dx);

duvdy = compactdiff(uv,dy,2);
% duvdy(:,1    ,:) = ( -uv(:,3  ,:) +4*uv(:, 2 ,:) -3*uv(:,1    ,:)) / (2*dy);
% duvdy(:,2:N-1,:) = (  uv(:,3:N,:)                -  uv(:,1:N-2,:)) / (2*dy);
% duvdy(:,    N,:) = (3*uv(:,  N,:) -4*uv(:,N-1,:) +  uv(:,  N-2,:)) / (2*dy);

dvvdy = compactdiff(vv,dy,2);
% dvvdy(:,1    ,:) = ( -vv(:,3  ,:) +4*vv(:, 2 ,:) -3*vv(:,1    ,:)) / (2*dy);
% dvvdy(:,2:N-1,:) = (  vv(:,3:N,:)                -  vv(:,1:N-2,:)) / (2*dy);
% dvvdy(:,    N,:) = (3*vv(:,  N,:) -4*vv(:,N-1,:) +  vv(:,  N-2,:)) / (2*dy);

% d2uudx2 = compactdiff2(uu,dx,1);
% % d2uudx2(1    ,:,:) = (uu(3  ,:,:) -2*uu(2    ,:,:) + uu(1    ,:,:)) / (dx^2);
% % d2uudx2(2:M-1,:,:) = (uu(3:M,:,:) -2*uu(2:M-1,:,:) + uu(1:M-2,:,:)) / (dx^2);
% % d2uudx2(M    ,:,:) = (uu(  M,:,:) -2*uu(  M-1,:,:) + uu(  M-2,:,:)) / (dx^2);
% 
% d2vvdy2 = compactdiff2(vv,dy,2);
% % d2vvdy2(:,1    ,:) = (vv(:,3  ,:) -2*vv(:,2    ,:) + vv(:,1    ,:)) / (dy^2);
% % d2vvdy2(:,2:N-1,:) = (vv(:,3:N,:) -2*vv(:,2:N-1,:) + vv(:,1:N-2,:)) / (dy^2);
% % d2vvdy2(:,    N,:) = (vv(:,  N,:) -2*vv(:,  N-1,:) + vv(:,  N-2,:)) / (dy^2);


%simpler to recalculate after BC's are done - use updated dudx and dvdy
% d2udxdt(:,:,1    ) = ( -dudx(:,:,3  ) +4*dudx(:,:, 2 ) -3*dudx(:,:,1    )) / (2*dt);
% d2udxdt(:,:,2:T-1) = (  dudx(:,:,3:T)                  -  dudx(:,:,1:T-2)) / (2*dt);
% d2udxdt(:,:,    T) = (3*dudx(:,:,  T) -4*dudx(:,:,T-1) +  dudx(:,:,  T-2)) / (2*dt);
% 
% d2vdydt(:,:,1    ) = ( -dvdy(:,:,3  ) +4*dvdy(:,:, 2 ) -3*dvdy(:,:,1    )) / (2*dt);
% d2vdydt(:,:,2:T-1) = (  dvdy(:,:,3:T)                  -  dvdy(:,:,1:T-2)) / (2*dt);
% d2vdydt(:,:,    T) = (3*dvdy(:,:,  T) -4*dvdy(:,:,T-1) +  dvdy(:,:,  T-2)) / (2*dt);

d2udydx = compactdiff(dudy,dx,1);
% d2udydx(1    ,:,:) = ( -dudy(3  ,:,:) +4*dudy( 2 ,:,:) -3*dudy(1    ,:,:)) / (2*dx);
% d2udydx(2:M-1,:,:) = (  dudy(3:M,:,:)                  -  dudy(1:M-2,:,:)) / (2*dx);
% d2udydx(M    ,:,:) = (3*dudy(M  ,:,:) -4*dudy(M-1,:,:) +  dudy(  M-2,:,:)) / (2*dx);

d2vdxdy = compactdiff(dvdx,dy,2);
% d2vdxdy(:,1    ,:) = ( -dvdx(:,3  ,:) +4*dvdx(:, 2 ,:) -3*dvdx(:,1    ,:)) / (2*dy);
% d2vdxdy(:,2:N-1,:) = (  dvdx(:,3:N,:)                  -  dvdx(:,1:N-2,:)) / (2*dy);
% d2vdxdy(:,    N,:) = (3*dvdx(:,  N,:) -4*dvdx(:,N-1,:) +  dvdx(:,  N-2,:)) / (2*dy);

% %calculate wall direction, normals, and derivatives at BC pts
% %still going to have problems in 1-node acute points
% %will crash if narrow channels close to domain wall face outwards (ie
% %must be more than 2 points existing in direction BC points towards)
% for n=1:length(bclist)
%     flat = 1e-6;
%     %disp(' xi  yi       s1       s2    nhat1    nhat2  N/S  E/W')
%     
%     %first point in list
%         xi=bclist{n}(1,1);
%         yi=bclist{n}(1,2);
%         s = bclist{n}(2,:) - bclist{n}(end,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         nmlist{n}(1,:) = nhat;  %store unit normal for plotting
%         %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
%         
%         if nhat(2) < -flat      % north
%             %use forward differences pointing downward
%             %fprintf('     N')
%             dudy  (xi, yi,:)    = (3*u(xi, yi,:) -4*u(xi, yi-1,:) + u(xi, yi-2,:)) / (2*dy);
%             dvdy  (xi, yi,:)    = (3*v(xi, yi,:) -4*v(xi, yi-1,:) + v(xi, yi-2,:)) / (2*dy);
%             d2udy2(xi, yi,:)    = (  u(xi, yi,:) -2*u(xi, yi-1,:) + u(xi, yi-2,:)) / (dy^2);
%             d2vdy2(xi, yi,:)    = (  v(xi, yi,:) -2*v(xi, yi-1,:) + v(xi, yi-2,:)) / (dy^2);
%             duvdy  (xi, yi,:)    = (3*uv(xi, yi,:) -4*uv(xi, yi-1,:) + uv(xi, yi-2,:)) / (2*dy);
%             dvvdy  (xi, yi,:)    = (3*vv(xi, yi,:) -4*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (2*dy);
%             d2vvdy2(xi, yi,:)    = (  vv(xi, yi,:) -2*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (dy^2);
%             
%         elseif nhat(2) > flat   % south
%             %use forward differences pointing upward
%             %fprintf('     S')
%             dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dy);
%             dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dy);
%             d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
%             d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
%             duvdy  (xi, yi,:) = (-uv(xi, yi+2,:) +4*uv(xi, yi+1,:) -3*uv(xi, yi,:)) / (2*dy);
%             dvvdy  (xi, yi,:) = (-vv(xi, yi+2,:) +4*vv(xi, yi+1,:) -3*vv(xi, yi,:)) / (2*dy);
%             d2vvdy2(xi, yi,:) = ( vv(xi, yi+2,:) -2*vv(xi, yi+1,:) +  vv(xi, yi,:)) / (dy^2);
% 
%         else                    % flat
%             %do nothing to y-derivatives
%             %use central difference centered on this point
%             %fprintf('     0')
%         end
% 
%         if nhat(1) > flat       %west
%             %use forward differences pointing right
%             %fprintf('     W\n')
%             dudx  (xi ,yi,:) = (-u(xi+2, yi,:) +4*u(xi+1, yi,:) -3*u(xi, yi,:)) / (2*dx);
%             dvdx  (xi ,yi,:) = (-v(xi+2, yi,:) +4*v(xi+1, yi,:) -3*v(xi, yi,:)) / (2*dx);
%             d2udx2(xi ,yi,:) = ( u(xi+2, yi,:) -2*u(xi+1, yi,:) +  u(xi, yi,:)) / (dx^2);
%             d2vdx2(xi ,yi,:) = ( v(xi+2, yi,:) -2*v(xi+1, yi,:) +  v(xi, yi,:)) / (dx^2);
%             duudx  (xi ,yi,:) = (-uu(xi+2, yi,:) +4*uu(xi+1, yi,:) -3*uu(xi, yi,:)) / (2*dx);
%             duvdx  (xi ,yi,:) = (-uv(xi+2, yi,:) +4*uv(xi+1, yi,:) -3*uv(xi, yi,:)) / (2*dx);
%             d2uudx2(xi ,yi,:) = ( uu(xi+2, yi,:) -2*uu(xi+1, yi,:) +  uu(xi, yi,:)) / (dx^2);
% 
%         elseif nhat(1) < -flat  %east
%             %use forward differences pointing left
%             %fprintf('     E\n')
%             dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
%             dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
%             d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
%             d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
%             duudx  (xi, yi,:) = (3*uu(xi, yi,:) -4*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (2*dx);
%             duvdx  (xi, yi,:) = (3*uv(xi, yi,:) -4*uv(xi-1, yi,:) + uv(xi-2, yi,:)) / (2*dx);
%             d2uudx2(xi, yi,:) = (  uu(xi, yi,:) -2*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (dx^2);
%             
%         else            
%             %do nothing to x-derivatives
%             %use central difference centered on this point
%             %fprintf('     0\n')
%         end
%         %pause
%         
%     %middle pts
%     for i=2:size(bclist{n},1)-1
%         xi=bclist{n}(i,1);
%         yi=bclist{n}(i,2);
%         s = bclist{n}(i+1,:) - bclist{n}(i-1,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         nmlist{n}(i,:) = nhat;  %store unit normal for plotting
%         %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
%         
%         if nhat(2) < -flat      % north
%             %use forward differences pointing downward
%             %fprintf('     N')
%             dudy  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi, yi-1,:) + u(xi, yi-2,:)) / (2*dy);
%             dvdy  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi, yi-1,:) + v(xi, yi-2,:)) / (2*dy);
%             d2udy2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi, yi-1,:) + u(xi, yi-2,:)) / (dy^2);
%             d2vdy2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi, yi-1,:) + v(xi, yi-2,:)) / (dy^2);
%             duvdy  (xi, yi,:)    = (3*uv(xi, yi,:) -4*uv(xi, yi-1,:) + uv(xi, yi-2,:)) / (2*dy);
%             dvvdy  (xi, yi,:)    = (3*vv(xi, yi,:) -4*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (2*dy);
%             d2vvdy2(xi, yi,:)    = (  vv(xi, yi,:) -2*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (dy^2);
%             
%         elseif nhat(2) > flat   % south
%             %use forward differences pointing upward
%             %fprintf('     S')
%             dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dy);
%             dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dy);
%             d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
%             d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
%             duvdy  (xi, yi,:) = (-uv(xi, yi+2,:) +4*uv(xi, yi+1,:) -3*uv(xi, yi,:)) / (2*dy);
%             dvvdy  (xi, yi,:) = (-vv(xi, yi+2,:) +4*vv(xi, yi+1,:) -3*vv(xi, yi,:)) / (2*dy);
%             d2vvdy2(xi, yi,:) = ( vv(xi, yi+2,:) -2*vv(xi, yi+1,:) +  vv(xi, yi,:)) / (dy^2);
%         else                    % flat
%             %do nothing to y-derivatives
%             %use central difference centered on this point
%             %fprintf('     0')
%         end
% 
%         if nhat(1) > flat       %west
%             %use forward differences pointing right
%             %fprintf('     W\n')
%             dudx  (xi ,yi,:) = (-u(xi+2, yi,:) +4*u(xi+1, yi,:) -3*u(xi, yi,:)) / (2*dx);
%             dvdx  (xi ,yi,:) = (-v(xi+2, yi,:) +4*v(xi+1, yi,:) -3*v(xi, yi,:)) / (2*dx);
%             d2udx2(xi ,yi,:) = ( u(xi+2, yi,:) -2*u(xi+1, yi,:) +  u(xi, yi,:)) / (dx^2);
%             d2vdx2(xi ,yi,:) = ( v(xi+2, yi,:) -2*v(xi+1, yi,:) +  v(xi, yi,:)) / (dx^2);
%             duudx  (xi ,yi,:) = (-uu(xi+2, yi,:) +4*uu(xi+1, yi,:) -3*uu(xi, yi,:)) / (2*dx);
%             duvdx  (xi ,yi,:) = (-uv(xi+2, yi,:) +4*uv(xi+1, yi,:) -3*uv(xi, yi,:)) / (2*dx);
%             d2uudx2(xi ,yi,:) = ( uu(xi+2, yi,:) -2*uu(xi+1, yi,:) +  uu(xi, yi,:)) / (dx^2);
% 
%         elseif nhat(1) < -flat  %east
%             %use forward differences pointing left
%             %fprintf('     E\n')
%             dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
%             dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
%             d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
%             d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
%             duudx  (xi, yi,:) = (3*uu(xi, yi,:) -4*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (2*dx);
%             duvdx  (xi, yi,:) = (3*uv(xi, yi,:) -4*uv(xi-1, yi,:) + uv(xi-2, yi,:)) / (2*dx);
%             d2uudx2(xi, yi,:) = (  uu(xi, yi,:) -2*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (dx^2);
%         else            
%             %do nothing to x-derivatives
%             %use central difference centered on this point
%             %fprintf('     0\n')
%         end
%         %pause
%     end %middle points
%         
%     %last point in list
%         xi=bclist{n}(end,1);
%         yi=bclist{n}(end,2);
%         s = bclist{n}(1,:) - bclist{n}(end-1,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         nmlist{n}(end+1,:) = nhat;  %store unit normal for plotting
%         %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
%         
%         if nhat(2) < -flat      % north
%             %use forward differences pointing downward
%             %fprintf('     N')
%             dudy  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi, yi-1,:) + u(xi, yi-2,:)) / (2*dy);
%             dvdy  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi, yi-1,:) + v(xi, yi-2,:)) / (2*dy);
%             d2udy2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi, yi-1,:) + u(xi, yi-2,:)) / (dy^2);
%             d2vdy2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi, yi-1,:) + v(xi, yi-2,:)) / (dy^2);
%             duvdy  (xi, yi,:)    = (3*uv(xi, yi,:) -4*uv(xi, yi-1,:) + uv(xi, yi-2,:)) / (2*dy);
%             dvvdy  (xi, yi,:)    = (3*vv(xi, yi,:) -4*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (2*dy);
%             d2vvdy2(xi, yi,:)    = (  vv(xi, yi,:) -2*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (dy^2);
%            
%         elseif nhat(2) > flat   % south
%             %use forward differences pointing upward
%             %fprintf('     S')
%             dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dy);
%             dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dy);
%             d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
%             d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
%             duvdy  (xi, yi,:) = (-uv(xi, yi+2,:) +4*uv(xi, yi+1,:) -3*uv(xi, yi,:)) / (2*dy);
%             dvvdy  (xi, yi,:) = (-vv(xi, yi+2,:) +4*vv(xi, yi+1,:) -3*vv(xi, yi,:)) / (2*dy);
%             d2vvdy2(xi, yi,:) = ( vv(xi, yi+2,:) -2*vv(xi, yi+1,:) +  vv(xi, yi,:)) / (dy^2);
%         else                    % flat
%             %do nothing to y-derivatives
%             %use central difference centered on this point
%             %fprintf('     0')
%         end
% 
%         if nhat(1) > flat       %west
%             %use forward differences pointing right
%             %fprintf('     W\n')
%             dudx  (xi ,yi,:) = (-u(xi+2, yi,:) +4*u(xi+1, yi,:) -3*u(xi, yi,:)) / (2*dx);
%             dvdx  (xi ,yi,:) = (-v(xi+2, yi,:) +4*v(xi+1, yi,:) -3*v(xi, yi,:)) / (2*dx);
%             d2udx2(xi ,yi,:) = ( u(xi+2, yi,:) -2*u(xi+1, yi,:) +  u(xi, yi,:)) / (dx^2);
%             d2vdx2(xi ,yi,:) = ( v(xi+2, yi,:) -2*v(xi+1, yi,:) +  v(xi, yi,:)) / (dx^2);
%             duudx  (xi ,yi,:) = (-uu(xi+2, yi,:) +4*uu(xi+1, yi,:) -3*uu(xi, yi,:)) / (2*dx);
%             duvdx  (xi ,yi,:) = (-uv(xi+2, yi,:) +4*uv(xi+1, yi,:) -3*uv(xi, yi,:)) / (2*dx);
%             d2uudx2(xi ,yi,:) = ( uu(xi+2, yi,:) -2*uu(xi+1, yi,:) +  uu(xi, yi,:)) / (dx^2);
% 
%         elseif nhat(1) < -flat  %east
%             %use forward differences pointing left
%             %fprintf('     E\n')
%             dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
%             dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
%             d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
%             d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
%             duudx  (xi, yi,:) = (3*uu(xi, yi,:) -4*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (2*dx);
%             duvdx  (xi, yi,:) = (3*uv(xi, yi,:) -4*uv(xi-1, yi,:) + uv(xi-2, yi,:)) / (2*dx);
%             d2uudx2(xi, yi,:) = (  uu(xi, yi,:) -2*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (dx^2);
%         else            
%             %do nothing to x-derivatives
%             %use central difference centered on this point
%             %fprintf('     0\n')
%         end
%         %pause
% %     figure(2345)
% %     plot(bclist{n}(:,1),bclist{n}(:,2))
% %     hold on
% %     quiver(bclist{n}(:,1),bclist{n}(:,2),nmlist{n}(:,1),nmlist{n}(:,2))
% %     hold off
%     %keyboard
% end %end recalculation of derivatives on BC

%dwdz correction based on continuity
dwdz = - (dudx + dvdy);
udwdz = u.*dwdz;
vdwdz = v.*dwdz;


% %use new values for 1st derivatives to calculate higher time derivatives
% d2udxdt = compactdiff(dudx,dt,3);
% % d2udxdt(:,:,1    ) = ( -dudx(:,:,3  ) +4*dudx(:,:, 2 ) -3*dudx(:,:,1    )) / (2*dt);
% % d2udxdt(:,:,2:T-1) = (  dudx(:,:,3:T)                  -  dudx(:,:,1:T-2)) / (2*dt);
% % d2udxdt(:,:,    T) = (3*dudx(:,:,  T) -4*dudx(:,:,T-1) +  dudx(:,:,  T-2)) / (2*dt);
% 
% d2vdydt = compactdiff(dvdy,dt,3);
% % d2vdydt(:,:,1    ) = ( -dvdy(:,:,3  ) +4*dvdy(:,:, 2 ) -3*dvdy(:,:,1    )) / (2*dt);
% % d2vdydt(:,:,2:T-1) = (  dvdy(:,:,3:T)                  -  dvdy(:,:,1:T-2)) / (2*dt);
% % d2vdydt(:,:,    T) = (3*dvdy(:,:,  T) -4*dvdy(:,:,T-1) +  dvdy(:,:,  T-2)) / (2*dt);
% 
% d_udwdz_dx = compactdiff(udwdz,dx,1);
% % d_udwdz_dx(1    ,:,:) = ( -udwdz(3  ,:,:) +4*udwdz( 2 ,:,:) -3*udwdz(1    ,:,:)) / (2*dx);
% % d_udwdz_dx(2:M-1,:,:) = (  udwdz(3:M,:,:)                   -  udwdz(1:M-2,:,:)) / (2*dx);
% % d_udwdz_dx(M    ,:,:) = (3*udwdz(M  ,:,:) -4*udwdz(M-1,:,:) +  udwdz(  M-2,:,:)) / (2*dx);
% 
% d_vdwdz_dy = compactdiff(vdwdz,dy,2);
% % d_vdwdz_dy(:,1    ,:) = ( -vdwdz(:,3  ,:) +4*vdwdz(:, 2 ,:) -3*vdwdz(:,1    ,:)) / (2*dy);
% % d_vdwdz_dy(:,2:N-1,:) = (  vdwdz(:,3:N,:)                   -  vdwdz(:,1:N-2,:)) / (2*dy);
% % d_vdwdz_dy(:,    N,:) = (3*vdwdz(:,  N,:) -4*vdwdz(:,N-1,:) +  vdwdz(:,  N-2,:)) / (2*dy);
% 
% d2uvdydx = compactdiff(duvdy,dx,1);
% % d2uvdydx(1    ,:,:) = ( -duvdy(3  ,:,:) +4*duvdy( 2 ,:,:) -3*duvdy(1    ,:,:)) / (2*dx);
% % d2uvdydx(2:M-1,:,:) = (  duvdy(3:M,:,:)                   -  duvdy(1:M-2,:,:)) / (2*dx);
% % d2uvdydx(M    ,:,:) = (3*duvdy(M  ,:,:) -4*duvdy(M-1,:,:) +  duvdy(  M-2,:,:)) / (2*dx);
% 
% d2uvdxdy = compactdiff(duvdx,dy,2);
% % d2uvdxdy(:,1    ,:) = ( -duvdx(:,3  ,:) +4*duvdx(:, 2 ,:) -3*duvdx(:,1    ,:)) / (2*dy);
% % d2uvdxdy(:,2:N-1,:) = (  duvdx(:,3:N,:)                   -  duvdx(:,1:N-2,:)) / (2*dy);
% % d2uvdxdy(:,    N,:) = (3*duvdx(:,  N,:) -4*duvdx(:,N-1,:) +  duvdx(:,  N-2,:)) / (2*dy);

d_dwdz_dx = compactdiff(dwdz,dx,1);
% d_dwdz_dx(1    ,:,:) = ( -dwdz(3  ,:,:) +4*dwdz( 2 ,:,:) -3*dwdz(1    ,:,:)) / (2*dx);
% d_dwdz_dx(2:M-1,:,:) = (  dwdz(3:M,:,:)                   - dwdz(1:M-2,:,:)) / (2*dx);
% d_dwdz_dx(M    ,:,:) = (3*dwdz(M  ,:,:) -4*dwdz(M-1,:,:) +  dwdz(  M-2,:,:)) / (2*dx);

d_dwdz_dy = compactdiff(dwdz,dy,2);
% d_dwdz_dy(:,1    ,:) = ( -dwdz(:,3  ,:) +4*dwdz(:, 2 ,:) -3*dwdz(:,1    ,:)) / (2*dy);
% d_dwdz_dy(:,2:N-1,:) = (  dwdz(:,3:N,:)                  -  dwdz(:,1:N-2,:)) / (2*dy);
% d_dwdz_dy(:,    N,:) = (3*dwdz(:,  N,:) -4*dwdz(:,N-1,:) +  dwdz(:,  N-2,:)) / (2*dy);


% %do higher order nested derivatives based on above BC derivatives
% %calculate wall direction, normals, and derivatives at BC pts
% %still going to have problems in 1-node acute points
% %will crash if narrow channels close to domain wall face outwards (ie
% %must be more than 2 points existing in direction BC points towards)
% for n=1:length(bclist)
%     flat = 1e-6;
%     %disp(' xi  yi       s1       s2    nhat1    nhat2  N/S  E/W')
%     
%     %first point in list
%         xi=bclist{n}(1,1);
%         yi=bclist{n}(1,2);
%         s = bclist{n}(2,:) - bclist{n}(end,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         nmlist{n}(1,:) = nhat;  %store unit normal for plotting
%         %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
%         
%         if nhat(2) < -flat      % north
%             %use forward differences pointing downward
%             %fprintf('     N')
%             d2vdxdy   (xi, yi,:) = (3*  dvdx(xi, yi,:) -4*  dvdx(xi, yi-1,:) +   dvdx(xi, yi-2,:)) / (2*dy);
% %             d3vdx2dy  (xi, yi,:) = (3*d2vdx2(xi, yi,:) -4*d2vdx2(xi, yi-1,:) + d2vdx2(xi, yi-2,:)) / (2*dy);
% %             d3vdy2dy  (xi, yi,:) = (3*d2vdy2(xi, yi,:) -4*d2vdy2(xi, yi-1,:) + d2vdy2(xi, yi-2,:)) / (2*dy);
%             d_dwdz_dy (xi, yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi, yi-1,:) +   dwdz(xi, yi-2,:)) / (2*dy);
% %             d_vdwdz_dy(xi, yi,:) = (3* vdwdz(xi, yi,:) -4* vdwdz(xi, yi-1,:) +  vdwdz(xi, yi-2,:)) / (2*dy);
% %             d2uvdxdy  (xi, yi,:) = (3* duvdx(xi, yi,:) -4* duvdx(xi, yi-1,:) +  duvdx(xi, yi-2,:)) / (2*dy);
%             
%         elseif nhat(2) > flat   % south
%             %use forward differences pointing upward
%             %fprintf('     S')
%             d2vdxdy   (xi, yi,:) = (-  dvdx(xi, yi+2,:) +4*  dvdx(xi, yi+1,:) -3*  dvdx(xi, yi,:)) / (2*dy);
% %             d3vdx2dy  (xi, yi,:) = (-d2vdx2(xi, yi+2,:) +4*d2vdx2(xi, yi+1,:) -3*d2vdx2(xi, yi,:)) / (2*dy);
% %             d3vdy2dy  (xi, yi,:) = (-d2vdy2(xi, yi+2,:) +4*d2vdy2(xi, yi+1,:) -3*d2vdy2(xi, yi,:)) / (2*dy);
%             d_dwdz_dy (xi, yi,:) = (-  dwdz(xi, yi+2,:) +4*  dwdz(xi, yi+1,:) -3*  dwdz(xi, yi,:)) / (2*dy);
% %             d_vdwdz_dy(xi, yi,:) = (- vdwdz(xi, yi+2,:) +4* vdwdz(xi, yi+1,:) -3* vdwdz(xi, yi,:)) / (2*dy);
% %             d2uvdxdy  (xi, yi,:) = (- duvdx(xi, yi+2,:) +4* duvdx(xi, yi+1,:) -3* duvdx(xi, yi,:)) / (2*dy);
%             
%         else                    % flat
%             %do nothing to y-derivatives
%             %use central difference centered on this point
%             %fprintf('     0')
%         end
% 
%         if nhat(1) > flat       %west
%             %use forward differences pointing right
%             %fprintf('     W\n')
%             d2udydx   (xi ,yi,:) = (-  dudy(xi+2, yi,:) +4*  dudy(xi+1, yi,:) -3*  dudy(xi, yi,:)) / (2*dx); 
% %             d3udx2dx  (xi ,yi,:) = (-d2udx2(xi+2, yi,:) +4*d2udx2(xi+1, yi,:) -3*d2udx2(xi, yi,:)) / (2*dx);
% %             d3udy2dx  (xi ,yi,:) = (-d2udy2(xi+2, yi,:) +4*d2udy2(xi+1, yi,:) -3*d2udy2(xi, yi,:)) / (2*dx);
%             d_dwdz_dx (xi ,yi,:) = (-  dwdz(xi+2, yi,:) +4*  dwdz(xi+1, yi,:) -3*  dwdz(xi, yi,:)) / (2*dx); 
% %             d_udwdz_dx(xi ,yi,:) = (- udwdz(xi+2, yi,:) +4* udwdz(xi+1, yi,:) -3* udwdz(xi, yi,:)) / (2*dx); 
% %             d2uvdydx  (xi ,yi,:) = (- duvdy(xi+2, yi,:) +4* duvdy(xi+1, yi,:) -3* duvdy(xi, yi,:)) / (2*dx); 
%             
%         elseif nhat(1) < -flat  %east
%             %use forward differences pointing left
%             %fprintf('     E\n')
%             d2udydx   (xi, yi,:) = (3*  dudy(xi, yi,:) -4*  dudy(xi-1, yi,:) +   dudy(xi-2, yi,:)) / (2*dx);
% %             d3udx2dx  (xi, yi,:) = (3*d2udx2(xi, yi,:) -4*d2udx2(xi-1, yi,:) + d2udx2(xi-2, yi,:)) / (2*dx);
% %             d3udy2dx  (xi, yi,:) = (3*d2udy2(xi, yi,:) -4*d2udy2(xi-1, yi,:) + d2udy2(xi-2, yi,:)) / (2*dx);
%             d_dwdz_dx (xi ,yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi-1, yi,:) +   dwdz(xi-2, yi,:)) / (2*dx); 
% %             d_udwdz_dx(xi ,yi,:) = (3* udwdz(xi, yi,:) -4* udwdz(xi-1, yi,:) +  udwdz(xi-2, yi,:)) / (2*dx); 
% %             d2uvdydx  (xi ,yi,:) = (3* duvdy(xi, yi,:) -4* duvdy(xi-1, yi,:) +  duvdy(xi-2, yi,:)) / (2*dx); 
%             
%         else            
%             %do nothing to x-derivatives
%             %use central difference centered on this point
%             %fprintf('     0\n')
%         end
%         %pause
%         
%     %middle pts
%     for i=2:size(bclist{n},1)-1
%         xi=bclist{n}(i,1);
%         yi=bclist{n}(i,2);
%         s = bclist{n}(i+1,:) - bclist{n}(i-1,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         nmlist{n}(i,:) = nhat;  %store unit normal for plotting
%         %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
%         
%         if nhat(2) < -flat      % north
%             %use forward differences pointing downward
%             %fprintf('     N')
%             d2vdxdy   (xi, yi,:) = (3*  dvdx(xi, yi,:) -4*  dvdx(xi, yi-1,:) +   dvdx(xi, yi-2,:)) / (2*dy);
% %             d3vdx2dy  (xi, yi,:) = (3*d2vdx2(xi, yi,:) -4*d2vdx2(xi, yi-1,:) + d2vdx2(xi, yi-2,:)) / (2*dy);
% %             d3vdy2dy  (xi, yi,:) = (3*d2vdy2(xi, yi,:) -4*d2vdy2(xi, yi-1,:) + d2vdy2(xi, yi-2,:)) / (2*dy);
%             d_dwdz_dy (xi, yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi, yi-1,:) +   dwdz(xi, yi-2,:)) / (2*dy);
% %             d_vdwdz_dy(xi, yi,:) = (3* vdwdz(xi, yi,:) -4* vdwdz(xi, yi-1,:) +  vdwdz(xi, yi-2,:)) / (2*dy);
% %             d2uvdxdy  (xi, yi,:) = (3* duvdx(xi, yi,:) -4* duvdx(xi, yi-1,:) +  duvdx(xi, yi-2,:)) / (2*dy);
%             
%         elseif nhat(2) > flat   % south
%             %use forward differences pointing upward
%             %fprintf('     S')
%             d2vdxdy   (xi, yi,:) = (-  dvdx(xi, yi+2,:) +4*  dvdx(xi, yi+1,:) -3*  dvdx(xi, yi,:)) / (2*dy);
% %             d3vdx2dy  (xi, yi,:) = (-d2vdx2(xi, yi+2,:) +4*d2vdx2(xi, yi+1,:) -3*d2vdx2(xi, yi,:)) / (2*dy);
% %             d3vdy2dy  (xi, yi,:) = (-d2vdy2(xi, yi+2,:) +4*d2vdy2(xi, yi+1,:) -3*d2vdy2(xi, yi,:)) / (2*dy);
%             d_dwdz_dy (xi, yi,:) = (-  dwdz(xi, yi+2,:) +4*  dwdz(xi, yi+1,:) -3*  dwdz(xi, yi,:)) / (2*dy);
% %             d_vdwdz_dy(xi, yi,:) = (- vdwdz(xi, yi+2,:) +4* vdwdz(xi, yi+1,:) -3* vdwdz(xi, yi,:)) / (2*dy);
% %             d2uvdxdy  (xi, yi,:) = (- duvdx(xi, yi+2,:) +4* duvdx(xi, yi+1,:) -3* duvdx(xi, yi,:)) / (2*dy);
%            
%         else                    % flat
%             %do nothing to y-derivatives
%             %use central difference centered on this point
%             %fprintf('     0')
%         end
% 
%         if nhat(1) > flat       %west
%             %use forward differences pointing right
%             %fprintf('     W\n')
%             d2udydx   (xi ,yi,:) = (-  dudy(xi+2, yi,:) +4*  dudy(xi+1, yi,:) -3*  dudy(xi, yi,:)) / (2*dx); 
% %             d3udx2dx  (xi ,yi,:) = (-d2udx2(xi+2, yi,:) +4*d2udx2(xi+1, yi,:) -3*d2udx2(xi, yi,:)) / (2*dx);
% %             d3udy2dx  (xi ,yi,:) = (-d2udy2(xi+2, yi,:) +4*d2udy2(xi+1, yi,:) -3*d2udy2(xi, yi,:)) / (2*dx);
%             d_dwdz_dx (xi ,yi,:) = (-  dwdz(xi+2, yi,:) +4*  dwdz(xi+1, yi,:) -3*  dwdz(xi, yi,:)) / (2*dx); 
% %             d_udwdz_dx(xi ,yi,:) = (- udwdz(xi+2, yi,:) +4* udwdz(xi+1, yi,:) -3* udwdz(xi, yi,:)) / (2*dx); 
% %             d2uvdydx  (xi ,yi,:) = (- duvdy(xi+2, yi,:) +4* duvdy(xi+1, yi,:) -3* duvdy(xi, yi,:)) / (2*dx); 
%            
%         elseif nhat(1) < -flat  %east
%             %use forward differences pointing left
%             %fprintf('     E\n')
%             d2udydx   (xi, yi,:) = (3*  dudy(xi, yi,:) -4*  dudy(xi-1, yi,:) +   dudy(xi-2, yi,:)) / (2*dx);
% %             d3udx2dx  (xi, yi,:) = (3*d2udx2(xi, yi,:) -4*d2udx2(xi-1, yi,:) + d2udx2(xi-2, yi,:)) / (2*dx);
% %             d3udy2dx  (xi, yi,:) = (3*d2udy2(xi, yi,:) -4*d2udy2(xi-1, yi,:) + d2udy2(xi-2, yi,:)) / (2*dx);
%             d_dwdz_dx (xi ,yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi-1, yi,:) +   dwdz(xi-2, yi,:)) / (2*dx); 
% %             d_udwdz_dx(xi ,yi,:) = (3* udwdz(xi, yi,:) -4* udwdz(xi-1, yi,:) +  udwdz(xi-2, yi,:)) / (2*dx); 
% %             d2uvdydx  (xi ,yi,:) = (3* duvdy(xi, yi,:) -4* duvdy(xi-1, yi,:) +  duvdy(xi-2, yi,:)) / (2*dx); 
%            
%         else            
%             %do nothing to x-derivatives
%             %use central difference centered on this point
%             %fprintf('     0\n')
%         end
%         %pause
%     end %middle points
%         
%     %last point in list
%         xi=bclist{n}(end,1);
%         yi=bclist{n}(end,2);
%         s = bclist{n}(1,:) - bclist{n}(end-1,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         nmlist{n}(end+1,:) = nhat;  %store unit normal for plotting
%         %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
%         
%         if nhat(2) < -flat      % north
%             %use forward differences pointing downward
%             %fprintf('     N')
%             d2vdxdy   (xi, yi,:) = (3*  dvdx(xi, yi,:) -4*  dvdx(xi, yi-1,:) +   dvdx(xi, yi-2,:)) / (2*dy);
% %             d3vdx2dy  (xi, yi,:) = (3*d2vdx2(xi, yi,:) -4*d2vdx2(xi, yi-1,:) + d2vdx2(xi, yi-2,:)) / (2*dy);
% %             d3vdy2dy  (xi, yi,:) = (3*d2vdy2(xi, yi,:) -4*d2vdy2(xi, yi-1,:) + d2vdy2(xi, yi-2,:)) / (2*dy);
%             d_dwdz_dy (xi, yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi, yi-1,:) +   dwdz(xi, yi-2,:)) / (2*dy);
% %             d_vdwdz_dy(xi, yi,:) = (3* vdwdz(xi, yi,:) -4* vdwdz(xi, yi-1,:) +  vdwdz(xi, yi-2,:)) / (2*dy);
% %             d2uvdxdy  (xi, yi,:) = (3* duvdx(xi, yi,:) -4* duvdx(xi, yi-1,:) +  duvdx(xi, yi-2,:)) / (2*dy);
%            
%         elseif nhat(2) > flat   % south
%             %use forward differences pointing upward
%             %fprintf('     S')
%             d2vdxdy   (xi, yi,:) = (-  dvdx(xi, yi+2,:) +4*  dvdx(xi, yi+1,:) -3*  dvdx(xi, yi,:)) / (2*dy);
% %             d3vdx2dy  (xi, yi,:) = (-d2vdx2(xi, yi+2,:) +4*d2vdx2(xi, yi+1,:) -3*d2vdx2(xi, yi,:)) / (2*dy);
% %             d3vdy2dy  (xi, yi,:) = (-d2vdy2(xi, yi+2,:) +4*d2vdy2(xi, yi+1,:) -3*d2vdy2(xi, yi,:)) / (2*dy);
%             d_dwdz_dy (xi, yi,:) = (-  dwdz(xi, yi+2,:) +4*  dwdz(xi, yi+1,:) -3*  dwdz(xi, yi,:)) / (2*dy);
% %             d_vdwdz_dy(xi, yi,:) = (- vdwdz(xi, yi+2,:) +4* vdwdz(xi, yi+1,:) -3* vdwdz(xi, yi,:)) / (2*dy);
% %             d2uvdxdy  (xi, yi,:) = (- duvdx(xi, yi+2,:) +4* duvdx(xi, yi+1,:) -3* duvdx(xi, yi,:)) / (2*dy);
%            
%         else                    % flat
%             %do nothing to y-derivatives
%             %use central difference centered on this point
%             %fprintf('     0')
%         end
% 
%         if nhat(1) > flat       %west
%             %use forward differences pointing right
%             %fprintf('     W\n')
%             d2udydx   (xi ,yi,:) = (-  dudy(xi+2, yi,:) +4*  dudy(xi+1, yi,:) -3*  dudy(xi, yi,:)) / (2*dx); 
% %             d3udx2dx  (xi ,yi,:) = (-d2udx2(xi+2, yi,:) +4*d2udx2(xi+1, yi,:) -3*d2udx2(xi, yi,:)) / (2*dx);
% %             d3udy2dx  (xi ,yi,:) = (-d2udy2(xi+2, yi,:) +4*d2udy2(xi+1, yi,:) -3*d2udy2(xi, yi,:)) / (2*dx);
%             d_dwdz_dx (xi ,yi,:) = (-  dwdz(xi+2, yi,:) +4*  dwdz(xi+1, yi,:) -3*  dwdz(xi, yi,:)) / (2*dx); 
% %             d_udwdz_dx(xi ,yi,:) = (- udwdz(xi+2, yi,:) +4* udwdz(xi+1, yi,:) -3* udwdz(xi, yi,:)) / (2*dx); 
% %             d2uvdydx  (xi ,yi,:) = (- duvdy(xi+2, yi,:) +4* duvdy(xi+1, yi,:) -3* duvdy(xi, yi,:)) / (2*dx); 
%             
%         elseif nhat(1) < -flat  %east
%             %use forward differences pointing left
%             %fprintf('     E\n')
%             d2udydx   (xi, yi,:) = (3*  dudy(xi, yi,:) -4*  dudy(xi-1, yi,:) +   dudy(xi-2, yi,:)) / (2*dx);
% %             d3udx2dx  (xi, yi,:) = (3*d2udx2(xi, yi,:) -4*d2udx2(xi-1, yi,:) + d2udx2(xi-2, yi,:)) / (2*dx);
% %             d3udy2dx  (xi, yi,:) = (3*d2udy2(xi, yi,:) -4*d2udy2(xi-1, yi,:) + d2udy2(xi-2, yi,:)) / (2*dx);
%             d_dwdz_dx (xi ,yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi-1, yi,:) +   dwdz(xi-2, yi,:)) / (2*dx); 
% %             d_udwdz_dx(xi ,yi,:) = (3* udwdz(xi, yi,:) -4* udwdz(xi-1, yi,:) +  udwdz(xi-2, yi,:)) / (2*dx); 
% %             d2uvdydx  (xi ,yi,:) = (3* duvdy(xi, yi,:) -4* duvdy(xi-1, yi,:) +  duvdy(xi-2, yi,:)) / (2*dx); 
%             
%         else            
%             %do nothing to x-derivatives
%             %use central difference centered on this point
%             %fprintf('     0\n')
%         end
%         %pause
% end %end recalculation of higher derivatives on BC

%this section was replaced - no longer check for ideal vs viscous
% %Calculate N-S field function        
% if strcmp(flowtype,'ideal')
%     f = - ( dudt + u.*dudx + v.*dudy );
%     g = - ( dvdt + u.*dvdx + v.*dvdy );
% elseif strcmp(flowtype,'viscous')
%     f = 1/Re * (d2udx2+d2udy2) - ( dudt + u.*dudx + v.*dudy );
%     g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + u.*dvdx + v.*dvdy );
% else
%     error('unknown flowtype');
% end
% % %Calculate N-S field function in conservative form      
% % if strcmp(flowtype,'ideal')
% %     f = - ( dudt + duudx + duvdy + udwdz);
% %     g = - ( dvdt + duvdx + dvvdy + vdwdz);
% % elseif strcmp(flowtype,'viscous')
% %     f = 1/Re * (d2udx2+d2udy2) - ( dudt + duudx + duvdy + udwdz);
% %     g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + duvdx + dvvdy + vdwdz);
% % else
% %     error('unknown flowtype');
% % end

%keyboard

%Calculate N-S field function        
if strcmp(solver,'conservative')
    f = 1/Re * (d2udx2+d2udy2) - ( dudt + duudx + duvdy + udwdz);
    g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + duvdx + dvvdy + vdwdz);
elseif strcmp(solver,'conservative2')
    f = 1/Re * (2*d2udx2 +   d2udy2 + d2vdxdy + d_dwdz_dx) - ( dudt + duudx + duvdy + udwdz);
    g = 1/Re * (  d2vdx2 + 2*d2vdy2 + d2udydx + d_dwdz_dy) - ( dvdt + duvdx + dvvdy + vdwdz);
elseif strcmp(solver,'standard')
    f = 1/Re * (d2udx2+d2udy2) - ( dudt + u.*dudx + v.*dudy );
    g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + u.*dvdx + v.*dvdy );
else
    error('unknown flowtype');
end





% %%%%%%%%%%
% %Calculate source terms
% %  should use approach from above on BC for this, or calculate each term
% %%%%%%%%%%
% dfdx(1    ,:,:) = (f(2  ,:,:) - f(1    ,:,:)) /    dx ;
% dfdx(2:M-1,:,:) = (f(3:M,:,:) - f(1:M-2,:,:)) / (2*dx);
% dfdx(M    ,:,:) = (f(M  ,:,:) - f(  M-1,:,:)) /    dx ;
% 
% dgdy(:,1    ,:) = (g(:,2  ,:) - g(:,1    ,:)) /    dy ;
% dgdy(:,2:N-1,:) = (g(:,3:N,:) - g(:,1:N-2,:)) / (2*dy);
% dgdy(:,    N,:) = (g(:,  N,:) - g(:,  N-1,:)) /    dy ;
% 
% Source = dfdx + dgdy;


%%Don't need source terms for line integrals
% % Source = -( dudx.*dudx + dudy.*dvdx ...
% %            +dvdx.*dudy + dvdy.*dvdy ); 
% 
% %assume continuity holds
% Source_old = -( dudx.*dudx + dudy.*dvdx ...
%              +dvdx.*dudy + dvdy.*dvdy ...
%              +dwdz.*dwdz );   
%        
% %use full form of RHS
% Source_new = 1/Re*( d3udx2dx + d3udy2dx + d3vdx2dy + d3vdy2dy) ...
%             - ( d2udxdt + d2vdydt ...
%                +dudx.*dudx + dudy.*dvdx ...
%                +dvdx.*dudy + dvdy.*dvdy ...
%                +u.*d2udx2 + v.*d2udydx + u.*d2vdxdy + v.*d2vdy2 );
%            
% %use conservative form of RHS
% Source = 1/Re*( d3udx2dx + d3udy2dx + d3vdx2dy + d3vdy2dy) ...
%             - ( d2udxdt + d2vdydt ...
%                +d2uudx2  + d2uvdydx + d_udwdz_dx ...
%                +d2uvdxdy + d2vvdy2  + d_vdwdz_dy );
%            
% %keyboard

%can probably prune the list of derivatives for line integrals - many terms
%are for the source term only
clear dudt dvdt dudx dvdx dudy dvdy d2udx2 d2udy2 d2vdx2 d2vdy2
clear dfdx dgdy uu uv vv udwdz vdwdz
clear d2udydx d3udx2dx d3udy2dx d2vdxdy d3vdx2dy d3vdy2dy d2udxdt d2vdydt
clear d2uudx2 d2uvdydx d_udwdz_dx d2uvdxdy d2vvdy2 d_vdwdz_dy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Don't need the coefficient matrix for line integrals
% % use SOC for middle of region
% % the problem has been discritized so that
% % Ap*phi_i,j    + 
% %  Ae*phi_i+1,j + Aw*phi_i-1,j  + 
% %  An*phi_i,j   + As*phi_i,j-1   =  S_i,j
% % 
% 
% %default discretization
% Ap(1:NI,1:NJ)   = -2*(1/dx^2+1/dy^2);
% Ae(1:NI,1:NJ)   = 1/dx^2;
% Aw(1:NI,1:NJ)   = 1/dx^2;
% An(1:NI,1:NJ)   = 1/dy^2;
% As(1:NI,1:NJ)   = 1/dy^2;
%  S(1:NI,1:NJ,:) = Source(1:NI,1:NJ,:);
% 
% 
% %calculate wall direction, normals, and BC
% for n=1:length(bclist)
% %     disp(' xi  yi       s1       s2    nhat1    nhat2       Ap       Ae       Aw       An       As        f        g        S')
%     
%     %first point in list
%         xi=bclist{n}(1,1);
%         yi=bclist{n}(1,2);
%         s = bclist{n}(2,:) - bclist{n}(end,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         
%         if nhat(1) > 0          %west
%             if nhat(2) > 0      % south
%                 Ap(xi, yi)   = -1/dx*nhat(1) - 1/dy*nhat(2);
%                 Ae(xi, yi)   =  1/dx*nhat(1);
%                 Aw(xi, yi)   =     0;
%                 An(xi, yi)   =                 1/dy*nhat(2);
%                 As(xi, yi)   =     0;
%             else                % north
%                 Ap(xi, yi)   = -1/dx*nhat(1) + 1/dy*nhat(2);
%                 Ae(xi, yi)   =  1/dx*nhat(1);
%                 Aw(xi, yi)   =     0;
%                 An(xi, yi)   =     0;
%                 As(xi, yi)   =               - 1/dy*nhat(2);   
%             end
%             
%         else                    %east
%             if nhat(2) > 0      % south
%                 Ap(xi, yi)   =  1/dx*nhat(1) - 1/dy*nhat(2);
%                 Ae(xi, yi)   =     0;
%                 Aw(xi, yi)   = -1/dx*nhat(1);
%                 An(xi, yi)   =                 1/dy*nhat(2);
%                 As(xi, yi)   =     0;
%             else                % north
%                 Ap(xi, yi)   =  1/dx*nhat(1) + 1/dy*nhat(2);
%                 Ae(xi, yi)   =     0;
%                 Aw(xi, yi)   = -1/dx*nhat(1);
%                 An(xi, yi)   =     0;
%                 As(xi, yi)   =                -1/dy*nhat(2);   
%             end
%         end    
%         S(xi, yi,:)  =  f(xi, yi,:)*nhat(1) + g(xi,yi,:)*nhat(2);
% 
% %         fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',xi, yi, s, nhat, Ap(xi,yi), Ae(xi,yi) ,Aw(xi,yi) ,An(xi,yi) ,As(xi,yi), f(xi,yi,1),  g(xi,yi,1), S(xi,yi,1))
% %         pause
%         
%     %middle pts
%     for i=2:size(bclist{n},1)-1
%         xi=bclist{n}(i,1);
%         yi=bclist{n}(i,2);
%         s = bclist{n}(i+1,:) - bclist{n}(i-1,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         
%         if nhat(1) > 0          %west
%             if nhat(2) > 0      % south
%                 Ap(xi, yi)   = -1/dx*nhat(1) - 1/dy*nhat(2);
%                 Ae(xi, yi)   =  1/dx*nhat(1);
%                 Aw(xi, yi)   =     0;
%                 An(xi, yi)   =                 1/dy*nhat(2);
%                 As(xi, yi)   =     0;
% 
%             else                % north
%                 Ap(xi, yi)   = -1/dx*nhat(1) + 1/dy*nhat(2);
%                 Ae(xi, yi)   =  1/dx*nhat(1);
%                 Aw(xi, yi)   =     0;
%                 An(xi, yi)   =     0;
%                 As(xi, yi)   =               - 1/dy*nhat(2);   
%             end
%             
%         else                    %east
%             if nhat(2) > 0      % south
%                 Ap(xi, yi)   =  1/dx*nhat(1) - 1/dy*nhat(2);
%                 Ae(xi, yi)   =     0;
%                 Aw(xi, yi)   = -1/dx*nhat(1);
%                 An(xi, yi)   =                 1/dy*nhat(2);
%                 As(xi, yi)   =     0;
%             else                % north
%                 Ap(xi, yi)   =  1/dx*nhat(1) + 1/dy*nhat(2);
%                 Ae(xi, yi)   =     0;
%                 Aw(xi, yi)   = -1/dx*nhat(1);
%                 An(xi, yi)   =     0;
%                 As(xi, yi)   =                -1/dy*nhat(2);   
%             end
%         end
%         S(xi, yi,:)  =  f(xi, yi,:)*nhat(1) + g(xi,yi,:)*nhat(2);
% %         fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',xi, yi, s, nhat, Ap(xi,yi), Ae(xi,yi) ,Aw(xi,yi) ,An(xi,yi) ,As(xi,yi), f(xi,yi,1),  g(xi,yi,1), S(xi,yi,1))
% %         pause        
%     end
% 
%     %last point in list
%         xi=bclist{n}(end,1);
%         yi=bclist{n}(end,2);
%         s = bclist{n}(1,:) - bclist{n}(end-1,:);
%         s = s/sqrt(s(1).^2+s(2).^2);
%         nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
%         
%         if nhat(1) > 0          %west
%             if nhat(2) > 0      % south
%                 Ap(xi, yi)   = -1/dx*nhat(1) - 1/dy*nhat(2);
%                 Ae(xi, yi)   =  1/dx*nhat(1);
%                 Aw(xi, yi)   =     0;
%                 An(xi, yi)   =                 1/dy*nhat(2);
%                 As(xi, yi)   =     0;
% 
%             else                % north
%                 Ap(xi, yi)   = -1/dx*nhat(1) + 1/dy*nhat(2);
%                 Ae(xi, yi)   =  1/dx*nhat(1);
%                 Aw(xi, yi)   =     0;
%                 An(xi, yi)   =     0;
%                 As(xi, yi)   =               - 1/dy*nhat(2);   
%             end
%             
%         else                    %east
%             if nhat(2) > 0      % south
%                 Ap(xi, yi)   =  1/dx*nhat(1) - 1/dy*nhat(2);
%                 Ae(xi, yi)   =     0;
%                 Aw(xi, yi)   = -1/dx*nhat(1);
%                 An(xi, yi)   =                 1/dy*nhat(2);
%                 As(xi, yi)   =     0;
%             else                % north
%                 Ap(xi, yi)   =  1/dx*nhat(1) + 1/dy*nhat(2);
%                 Ae(xi, yi)   =     0;
%                 Aw(xi, yi)   = -1/dx*nhat(1);
%                 An(xi, yi)   =     0;
%                 As(xi, yi)   =                -1/dy*nhat(2);   
%             end
%         end
%         S(xi, yi,:)  =  f(xi, yi,:)*nhat(1) + g(xi,yi,:)*nhat(2);
% %         fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',xi, yi, s, nhat, Ap(xi,yi), Ae(xi,yi) ,Aw(xi,yi) ,An(xi,yi) ,As(xi,yi), f(xi,yi,1),  g(xi,yi,1), S(xi,yi,1))
% %         pause
% end %for n=1:length(bclist) 
% 
% % %Bad idea - convergence is terrible - just cheat and shift the 
% % if nargin>=9
% %     %Fixed pressure measurement at a point
% %     Ap(Px,Py)   = 1;
% %     Ae(Px,Py)   = 0;
% %     Aw(Px,Py)   = 0;
% %     An(Px,Py)   = 0;
% %     As(Px,Py)   = 0;
% %      S(Px,Py,:) = Pref;
% % end
% 
% %convert coefficient lists to a sparse matrix
% AA = spdiags( [reshape(As,NI*NJ,1), ...
%                reshape(Aw,NI*NJ,1), ...
%                reshape(Ap,NI*NJ,1), ...
%                reshape(Ae,NI*NJ,1), ...
%                reshape(An,NI*NJ,1)] , [+NI,+1,0,-1,-NI], NI*NJ, NJ*NI ).';

% %need this?
% %solution copies previous solution to initial guess on each iteration, so
% %we need to make sure there is something in 2 and 3
% phi(:,:,2)=phi(:,:,3);  %fill in initial values for 

ktot = 0;
%figure(942),hold off,semilogy(0,0,'w')

% %garbage line integral tests
% temp1 = zeros(NI,NJ);
% temp2 = zeros(NI,NJ);
% temp3 = zeros(NI,NJ);
% temp4 = zeros(NI,NJ);
% for i=1:1:NI ,for j = 1:1:NJ ,temp1(i,j) = trapz(x( 1:1:i ,   1   ),f( 1:1:i ,    1  ,5),1) +trapz(y(   i   , 1:1:j ),g(   i   , 1:1:j ,5),2);end,end
% for i=NI:-1:1,for j = 1:1:NJ ,temp2(i,j) = trapz(y(  NI   , 1:1:j ),g(   NI  , 1:1:j ,5),2) +trapz(x(NI:-1:i,   j   ),f(NI:-1:i,   j   ,5),1);end,end
% for i=NI:-1:1,for j = NJ:-1:1,temp3(i,j) = trapz(x(NI:-1:i,   NJ  ),f(NI:-1:i,   NJ  ,5),1) +trapz(y(   i   ,NJ:-1:j),g(   i   ,NJ:-1:j,5),2);end,end
% for i=1:1:NI ,for j = NJ:-1:1,temp4(i,j) = trapz(y(   1   ,NJ:-1:j),g(   1   ,NJ:-1:j,5),2) +trapz(x( 1:1:i ,   j   ),f( 1:1:i ,   j   ,5),1);end,end
% 
% figure(801),surf(x,y,temp1)
% figure(802),surf(x,y,temp2)
% figure(803),surf(x,y,temp3)
% figure(804),surf(x,y,temp4)
% 
% tempavg = (temp1+temp2+temp3+temp4)/4;
% figure(805),surf(x,y,tempavg)

% keyboard
% %end garbage
% 
% return
pause(0.001)

% %startup condition to stop 1st timestep from converging
% fprintf('initializing first timestep... ')
%     [phi,flag,relres,iter,resvec] = gmres(AA, reshape(S(:,:,1),NI*NJ,1), restarts, tol, max_iterations, [], [], phi);
% fprintf('done\n')
% %fprintf('done, (%7.3f)\n',toc)

% fprintf('Performing incomplete LU factorization... ')
% [L,U]=luinc(AA,1e-5);
% toc

%solve for pressures
fprintf('     t       k        resid  flag        time\n')
fprintf('  ----  ------  -----------  ----  ----------\n')
% fprintf('  %4i  %6i  %11.4g  %4i  %10.3f\n',0,(iter(1)-1)*restarts+iter(2),relres,flag,toc)

%keyboard


for t=1:NT
    p1 = zeros(NI,NJ);
    p2 = zeros(NI,NJ);
    p3 = zeros(NI,NJ);
    p4 = zeros(NI,NJ);
    p5 = zeros(NI,NJ);
    p6 = zeros(NI,NJ);
    p7 = zeros(NI,NJ);
    p8 = zeros(NI,NJ);
    %dd = sqrt(dx.^2 + dy.^2);
    
    %from SW->SE corner
    j=1;
        for i=2:NI
            p1(i,j) =    p1(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx; %from W
        end
    for j=2:NJ
        i = 1;
            p1(i,j) = (  p1(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy ...  %from S
                       + p1(i+1,j-1) - 1/2*(f(i+1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i+1,j-1,t)+g(i,j,t))*dy )/2; %from SE
        for i=2:NI-1
            p1(i,j) = (  p1(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p1(i-1,j-1) + 1/2*(f(i-1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i-1,j-1,t)+g(i,j,t))*dy ...  %from SW
                       + p1(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy ...; %from S
                       + p1(i+1,j-1) - 1/2*(f(i+1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i+1,j-1,t)+g(i,j,t))*dy )/4; %from SE
        end
        i = NI;
            p1(i,j) = (  p1(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p1(i-1,j-1) + 1/2*(f(i-1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i-1,j-1,t)+g(i,j,t))*dy ...  %from SW
                       + p1(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy )/3; %from S
    end

    
    %from SW->NW corner
    i=1;
        for j=2:NJ
            p5(i,j) =    p5(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy;     %from S
        end
    for i=2:NI
        j = 1;
            p5(i,j) = (  p5(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p5(i-1,j+1) + 1/2*(f(i-1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i-1,j+1,t)+g(i,j,t))*dy )/2; %from NW
        for j=2:NJ-1
            p5(i,j) = (  p5(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p5(i-1,j-1) + 1/2*(f(i-1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i-1,j-1,t)+g(i,j,t))*dy ...  %from SW
                       + p5(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy ...; %from S
                       + p5(i-1,j+1) + 1/2*(f(i-1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i-1,j+1,t)+g(i,j,t))*dy )/4; %from NW
        end
        j = NJ;
            p5(i,j) = (  p5(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p5(i-1,j-1) + 1/2*(f(i-1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i-1,j-1,t)+g(i,j,t))*dy ...  %from SW
                       + p5(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy )/3; %from S
    end
    
    %from SE->SW corner
    j=1;
        for i=NI-1:-1:1
            p2(i,j) =    p2(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx; %from E
        end
    for j=2:NJ
        i = NI;
            p2(i,j) = (  p2(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy ...  %from S
                       + p2(i-1,j-1) - 1/2*(f(i-1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i-1,j-1,t)+g(i,j,t))*dy )/2; %from SW
        for i=NI-1:-1:2
            p2(i,j) = (  p2(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p2(i+1,j-1) - 1/2*(f(i+1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i+1,j-1,t)+g(i,j,t))*dy ...  %from SE
                       + p2(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy ...; %from S
                       + p2(i-1,j-1) + 1/2*(f(i-1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i-1,j-1,t)+g(i,j,t))*dy )/4; %from SW
        end
        i = 1;
            p2(i,j) = (  p2(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p2(i+1,j-1) - 1/2*(f(i+1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i+1,j-1,t)+g(i,j,t))*dy ...  %from SE
                       + p2(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy )/3; %from S
    end
    
    %from SE->NE corner
    i=NI;
        for j=2:1:NJ
            p6(i,j) =    p6(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy;     %from S
        end
    for i=NI-1:-1:1
        j = 1;
            p6(i,j) = (  p6(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p6(i+1,j+1) - 1/2*(f(i+1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i+1,j+1,t)+g(i,j,t))*dy )/2; %from NE
        for j=2:1:NJ-1
            p6(i,j) = (  p6(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p6(i+1,j+1) - 1/2*(f(i+1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i+1,j+1,t)+g(i,j,t))*dy ...  %from NE
                       + p6(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy ...  %from S
                       + p6(i+1,j-1) - 1/2*(f(i+1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i+1,j-1,t)+g(i,j,t))*dy )/4; %from SE
        end
        j = NJ;
            p6(i,j) = (  p6(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p6(i+1,j-1) - 1/2*(f(i+1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i+1,j-1,t)+g(i,j,t))*dy ...  %from SE
                       + p6(i  ,j-1)                                  + 1/2*(g(i  ,j-1,t)+g(i,j,t))*dy )/3; %from S
    end

    %from NE->NW corner
    j=NJ;
        for i=NI-1:-1:1
            p3(i,j) =    p3(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx;                                      %from E
        end
    for j=NJ-1:-1:1
        i = NI;
            p3(i,j) = (  p3(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy ...  %from N
                       + p3(i-1,j+1) + 1/2*(f(i-1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i-1,j+1,t)+g(i,j,t))*dy )/2; %from NW
        for i=NI-1:-1:2
            p3(i,j) = (  p3(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p3(i+1,j+1) - 1/2*(f(i+1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i+1,j+1,t)+g(i,j,t))*dy ...  %from NE
                       + p3(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy ...  %from N
                       + p3(i-1,j+1) + 1/2*(f(i-1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i-1,j+1,t)+g(i,j,t))*dy )/4; %from NW
        end
        i = 1;
            p3(i,j) = (  p3(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p3(i+1,j+1) - 1/2*(f(i+1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i+1,j+1,t)+g(i,j,t))*dy ...  %from NE
                       + p3(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy )/3; %from N
    end

    %from NE->SE corner
    i=NI;
        for j=NJ-1:-1:1
            p7(i,j) =    p7(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy;     %from N
        end
    for i=NI-1:-1:1
        j = NJ;
            p7(i,j) = (  p7(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p7(i+1,j-1) - 1/2*(f(i+1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i+1,j-1,t)+g(i,j,t))*dy )/2; %from SE
        for j=NJ-1:-1:2
            p7(i,j) = (  p7(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p7(i+1,j+1) - 1/2*(f(i+1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i+1,j+1,t)+g(i,j,t))*dy ...  %from NE
                       + p7(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy ...  %from N
                       + p7(i+1,j-1) - 1/2*(f(i+1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i+1,j-1,t)+g(i,j,t))*dy )/4; %from SE
        end
        j = 2;
            p7(i,j) = (  p7(i+1,j  ) - 1/2*(f(i+1,j  ,t)+f(i,j,t))*dx                                  ...  %from E
                       + p7(i+1,j+1) - 1/2*(f(i+1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i+1,j+1,t)+g(i,j,t))*dy ...  %from NE
                       + p7(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy )/3; %from N
    end
    
    %from NW->SW corner
    i=1;
        for j=NJ-1:-1:1
            p4(i,j) =    p4(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy;     %from N
        end
    for i=2:1:NI
        j = NJ;
            p4(i,j) = (  p4(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p4(i-1,j-1) + 1/2*(f(i-1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i-1,j-1,t)+g(i,j,t))*dy )/2; %from SW
        for j=NJ-1:-1:2
            p4(i,j) = (  p4(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p4(i-1,j-1) + 1/2*(f(i-1,j-1,t)+f(i,j,t))*dx + 1/2*(g(i-1,j-1,t)+g(i,j,t))*dy ...  %from SW
                       + p4(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy ...  %from N
                       + p4(i-1,j+1) + 1/2*(f(i-1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i-1,j+1,t)+g(i,j,t))*dy )/4; %from NW
        end
        j = 2;
            p4(i,j) = (  p4(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p4(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy ...  %from N
                       + p4(i-1,j+1) + 1/2*(f(i-1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i-1,j+1,t)+g(i,j,t))*dy )/3; %from NW
    end   
    
    %from NW->NE corner
    j=NJ;
        for i=2:1:NI
            p8(i,j) =    p8(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx;                                      %from W
        end
    for j=NJ-1:-1:1
        i = 1;
            p8(i,j) = (  p8(i+1,j+1) - 1/2*(f(i+1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i+1,j+1,t)+g(i,j,t))*dy ...  %from NE
                       + p8(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy )/2; %from N
        for i=2:1:NI-1
            p8(i,j) = (  p8(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p8(i+1,j+1) - 1/2*(f(i+1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i+1,j+1,t)+g(i,j,t))*dy ...  %from NE
                       + p8(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy ...  %from N
                       + p8(i-1,j+1) + 1/2*(f(i-1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i-1,j+1,t)+g(i,j,t))*dy )/4; %from NW
        end
        i = NI;
            p8(i,j) = (  p8(i-1,j  ) + 1/2*(f(i-1,j  ,t)+f(i,j,t))*dx                                  ...  %from W
                       + p8(i  ,j+1)                                  - 1/2*(g(i  ,j+1,t)+g(i,j,t))*dy ...  %from N
                       + p8(i-1,j+1) + 1/2*(f(i-1,j+1,t)+f(i,j,t))*dx - 1/2*(g(i-1,j+1,t)+g(i,j,t))*dy )/3; %from NW
    end   
   
%     figure(801),surf(x,y,p1(:,:)),title('p1: SW->SE')%,pause(0.0001)
%     figure(802),surf(x,y,p2(:,:)),title('p2: SE->SW')%,pause(0.0001)
%     figure(803),surf(x,y,p3(:,:)),title('p3: NE->NW')%,pause(0.0001)
%     figure(804),surf(x,y,p3(:,:)),title('p4: NW->SW')%,pause(0.0001)
%     figure(805),surf(x,y,p5(:,:)),title('p5: SW->NW')%,pause(0.0001)
%     figure(806),surf(x,y,p6(:,:)),title('p6: SE->NE')%,pause(0.0001)
%     figure(807),surf(x,y,p7(:,:)),title('p7: NE->SE')%,pause(0.0001)
%     figure(808),surf(x,y,p8(:,:)),title('p8: NW->NE'),pause(0.0001)
    
    flag=0;relres=0;iter=[1,1];,resvec=0;
    %[phi,flag,relres,iter,resvec] = gmres(AA, reshape(S(:,:,t),NI*NJ,1), restarts, tol, max_iterations, L, U, phi);
    %[phi]=AA\reshape(S(:,:,t),NI*NJ,1);  flag=0;relres=0;iter=[1,1];,resvec=0;

    %store final solution into pressure table
    %p(:,:,t) = reshape(phi,NI,NJ);
    p(:,:,t) = (p1+p2+p3+p4+p5+p6+p7+p8)/8;
    
    %use an offset to peg the reference pressure over time... 
    %is this okay to do?
    if nargin>=9
        %pause
        deltaP = Pref(t) - p(Px,Py,t);
        p(:,:,t) = p(:,:,t) + deltaP * ones(NI,NJ);
    end

    %convergence information
    fprintf('  %4i  %6i  %11.4g  %4i  %10.3f\n',t,(iter(1)-1)*restarts+iter(2),relres,flag,toc)
%     %plot current result for testing
%     figure(944),surf(x,y,p(:,:,t)),pause(0.0001)
%     figure(945),contourf(x,y,p(:,:,t),20),pause(0.0001)
%     figure(944),mesh(x,y,p(:,:,t)),pause(0.0001),colorbar
%     figure(945),imagesc([x(1,1),x(end,end)],[y(1,1),y(end,end)],p(:,:,t)),pause(0.0001),axis xy, colorbar
%     %plot convergence
%     normS = norm(reshape(S(:,:,t),NI*NJ,1));
%     figure(943),semilogy(resvec/normS),title(['residual for timestep t=',num2str(t)])
%     figure(942),hold on,semilogy(ktot+[0:length(resvec)-1],resvec/normS),title(['residual for timestep t=',num2str(t)])
    pause(1e-4)
    
%     figure(942),hold on,semilogy(ktot+[2:k],rho_rms(2:k)),title(['rho_r_m_s for timestep t=',num2str(t)])
%     hold off
    %figure(943),plot(rho_rms),title(['e_r_m_s for timestep t=',num2str(t)])
    %pause
       
     ktot = ktot + length(resvec)-1;
    

end %proceed to next time step


% figure(1)
% semilogy(1:k,rho_rms)
% xlabel('iteration,k'),ylabel('rho_r_m_s'),hold on
% title('convergence')
% 
% figure(2)
% semilogy(1:k,e_rms)
% xlabel('iteration,k'),ylabel('e_r_m_s'),hold on
% 
% figure(2000+1/d)
% contourf(x,y,phi(:,:,2),(0:.1:1)), colorbar, xlabel('x'), ylabel('y')
% title(strcat('ADI solution for dx=dy=',num2str(d),', k=',num2str(k)))
% 
% figure(2100+1/d)
% resid = (abs(phi_exact - phi(:,:,2))./phi_exact)*100;
% contour(x,y,resid); colorbar,xlabel('x'), ylabel('y')
% title(strcat('% resid for ADI solution at dx=dy=',num2str(d),', k=',num2str(k)))
% 
% figure(2200+1/d)
% contourf(x,y,phi_exact,(0:.1:1))
% title('exact solution'), colorbar, xlabel('x'), ylabel('y')

% catch
%     disp(lasterr)
%     keyboard
% end
