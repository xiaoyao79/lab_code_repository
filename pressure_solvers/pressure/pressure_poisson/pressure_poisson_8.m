function [p,S]=pressure_poisson_8(x,y,u,v,Re,dt,solver,bclist,Pref,Px,Py)
%function [S]=pressure_solve(x,y,u,v,Re,dt,[solver],[bclist],[Pref, Px, Py])
%
% Returns p(x,y,t) from u(x,y,t) and v(x,y,t) using the pressure Poisson 
% equation along with Navier-Stokes equations for 2D viscous flow:
%   del^2(P) = S(x,y)
%     S(x,y) = df/dx + dg/dy
%          f = dP/dx
%          g = dP/dy
%
% SOLVER must be either "conservative", "standard", or "simple". The 
% default is a conservative formulation for the convective terms:
%   "standard" Form:
%     dP/dx = 1/Re * (d2u/dx2+d2u/dy2) - ( du/dt + u*du/dx + v*du/dy )
%     dP/dy = 1/Re * (d2v/dx2+d2v/dy2) - ( dv/dt + u*dv/dx + v*dv/dy )
%
%   "conservative" Form:
%     dP/dx = 1/Re * (d2v/dx2+d2v/dy2) - (du/dt+d/dx[uu]+d/dy[uv]+u*dw/dz)
%     dP/dy = 1/Re * (d2v/dx2+d2v/dy2) - (dv/dt+d/dx[vu]+d/dy[vv]+v*dw/dz)
%     dw/dz = -( du/dx + dv/dy )
%
%   but for incompressible flow, the source term simplifies to:
%   "simple"
%      S(x,y) = du_i/dx_j*du_j/dx_i
%             = (du/dx)^2 + 2*du/dy*dv/dx + (dv/dy)^2
%
% Assume that X and Y are uniformly spaced, but dX and dY do not have 
% to be the same.
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

% changed to compact scheme for source terms - breaks uneven boundaries
% updated derivative calc on boundaries of domains...
% updated derivative calc on given boundary conditions
% added pressure Dirichlet BC Pref in interior at Px and Py
% using direct backsolve to calculate the pressures - ADI was poor, and
%  GMRES with LU decomposition worked, but was slower
% using incomplete LU factorization is waaay faster, but just using AA\S is
%  faster still - so I switched
% removed ideal flow option, replaced with switch for 3 forms of source
%  term
% pressure_poisson_1 was based off an out-of-date code branch - this is
%  forked from pressure_lines_2
% carried correction on d2vdxdy from pressure_direct2.m
% changed Pref calculation back to offset, from direct forcing in source
% added check to see if enough time steps (>=7) for compactdiff, else go
%  back to central differences
% changed back to central differences for derivatives - reactivated BC's
% v7 - for reprocessing for journal version.  Add 1 dirichlet point to
%    preserve stability.

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
    solver = 'conservative';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
tic
restarts        = 20;
max_iterations  = 50;
tol             = 1e-5;

dx = x(1,2) - x(1,1);
dy = y(2,1) - y(1,1);

% storage for source term
S = zeros(NI,NJ,NT);

%discretization defined on the following neighboring points
AA = sparse(NI*NJ,NI*NJ);
An = zeros(NI,NJ);  
As = zeros(NI,NJ);
Aw = zeros(NI,NJ);
Ae = zeros(NI,NJ);
Ap = zeros(NI,NJ);

%store pressure in phi during solution
phi = zeros(NI*NJ,1);

%storage for final pressure solutions
p = zeros(NI,NJ,NT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the source term 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_size = size(u);   %u and v should be the same size
M = u_size(1);      %number of points in x
N = u_size(2);      %number of points in y
T = u_size(3);      %number of points in t

fprintf('finding basic derivatives... ')

uu = u.*u;
uv = u.*v;
vv = v.*v;

    dudt(:,:,1    ) = ( -u(:,:,3  ) +4*u(:,:, 2 ) -3*u(:,:,1    )) / (2*dt);
    dudt(:,:,2:T-1) = (  u(:,:,3:T)               -  u(:,:,1:T-2)) / (2*dt);
    dudt(:,:,    T) = (3*u(:,:,  T) -4*u(:,:,T-1) +  u(:,:,  T-2)) / (2*dt);
    
    dvdt(:,:,1    ) = ( -v(:,:,3  ) +4*v(:,:, 2 ) -3*v(:,:,1    )) / (2*dt);
    dvdt(:,:,2:T-1) = (  v(:,:,3:T)               -  v(:,:,1:T-2)) / (2*dt);
    dvdt(:,:,    T) = (3*v(:,:,  T) -4*v(:,:,T-1) +  v(:,:,  T-2)) / (2*dt);

dudx(1    ,:,:) = ( -u(3  ,:,:) +4*u( 2 ,:,:) -3*u(1    ,:,:)) / (2*dx);
dudx(2:M-1,:,:) = (  u(3:M,:,:)               -  u(1:M-2,:,:)) / (2*dx);
dudx(M    ,:,:) = (3*u(M  ,:,:) -4*u(M-1,:,:) +  u(  M-2,:,:)) / (2*dx);

dvdx(1    ,:,:) = ( -v(3  ,:,:) +4*v( 2 ,:,:) -3*v(1    ,:,:)) / (2*dx);
dvdx(2:M-1,:,:) = (  v(3:M,:,:)               -  v(1:M-2,:,:)) / (2*dx);
dvdx(M    ,:,:) = (3*v(M  ,:,:) -4*v(M-1,:,:) +  v(  M-2,:,:)) / (2*dx);

dudy(:,1    ,:) = ( -u(:,3  ,:) +4*u(:, 2 ,:) -3*u(:,1    ,:)) / (2*dy);
dudy(:,2:N-1,:) = (  u(:,3:N,:)               -  u(:,1:N-2,:)) / (2*dy);
dudy(:,    N,:) = (3*u(:,  N,:) -4*u(:,N-1,:) +  u(:,  N-2,:)) / (2*dy);

dvdy(:,1    ,:) = ( -v(:,3  ,:) +4*v(:, 2 ,:) -3*v(:,1    ,:)) / (2*dy);
dvdy(:,2:N-1,:) = (  v(:,3:N,:)               -  v(:,1:N-2,:)) / (2*dy);
dvdy(:,    N,:) = (3*v(:,  N,:) -4*v(:,N-1,:) +  v(:,  N-2,:)) / (2*dy);

d2udx2(1    ,:,:) = (u(3  ,:,:) -2*u(2    ,:,:) + u(1    ,:,:)) / (dx^2);
d2udx2(2:M-1,:,:) = (u(3:M,:,:) -2*u(2:M-1,:,:) + u(1:M-2,:,:)) / (dx^2);
d2udx2(M    ,:,:) = (u(  M,:,:) -2*u(  M-1,:,:) + u(  M-2,:,:)) / (dx^2);

d2udy2(:,1    ,:) = (u(:,3  ,:) -2*u(:,2    ,:) + u(:,1    ,:)) / (dy^2);
d2udy2(:,2:N-1,:) = (u(:,3:N,:) -2*u(:,2:N-1,:) + u(:,1:N-2,:)) / (dy^2);
d2udy2(:,    N,:) = (u(:,  N,:) -2*u(:,  N-1,:) + u(:,  N-2,:)) / (dy^2);

d2vdx2(1    ,:,:) = (v(3  ,:,:) -2*v(2    ,:,:) + v(1    ,:,:)) / (dx^2);
d2vdx2(2:M-1,:,:) = (v(3:M,:,:) -2*v(2:M-1,:,:) + v(1:M-2,:,:)) / (dx^2);
d2vdx2(M    ,:,:) = (v(M  ,:,:) -2*v(  M-1,:,:) + v(  M-2,:,:)) / (dx^2);

d2vdy2(:,1    ,:) = (v(:,3  ,:) -2*v(:,2    ,:) + v(:,1    ,:)) / (dy^2);
d2vdy2(:,2:N-1,:) = (v(:,3:N,:) -2*v(:,2:N-1,:) + v(:,1:N-2,:)) / (dy^2);
d2vdy2(:,    N,:) = (v(:,  N,:) -2*v(:,  N-1,:) + v(:,  N-2,:)) / (dy^2);

d3udx2dx(1    ,:,:) = ( -d2udx2(3  ,:,:) +4*d2udx2( 2 ,:,:) -3*d2udx2(1    ,:,:)) / (2*dx);
d3udx2dx(2:M-1,:,:) = (  d2udx2(3:M,:,:)                    -  d2udx2(1:M-2,:,:)) / (2*dx);
d3udx2dx(M    ,:,:) = (3*d2udx2(M  ,:,:) -4*d2udx2(M-1,:,:) +  d2udx2(  M-2,:,:)) / (2*dx);

d3udy2dx(1    ,:,:) = ( -d2udy2(3  ,:,:) +4*d2udy2( 2 ,:,:) -3*d2udy2(1    ,:,:)) / (2*dx);
d3udy2dx(2:M-1,:,:) = (  d2udy2(3:M,:,:)                    -  d2udy2(1:M-2,:,:)) / (2*dx);
d3udy2dx(M    ,:,:) = (3*d2udy2(M  ,:,:) -4*d2udy2(M-1,:,:) +  d2udy2(  M-2,:,:)) / (2*dx);

d3vdx2dy(:,1    ,:) = ( -d2vdx2(:,3  ,:) +4*d2vdx2(:, 2 ,:) -3*d2vdx2(:,1    ,:)) / (2*dy);
d3vdx2dy(:,2:N-1,:) = (  d2vdx2(:,3:N,:)                    -  d2vdx2(:,1:N-2,:)) / (2*dy);
d3vdx2dy(:,    N,:) = (3*d2vdx2(:,  N,:) -4*d2vdx2(:,N-1,:) +  d2vdx2(:,  N-2,:)) / (2*dy);

d3vdy2dy(:,1    ,:) = ( -d2vdy2(:,3  ,:) +4*d2vdy2(:, 2 ,:) -3*d2vdy2(:,1    ,:)) / (2*dy);
d3vdy2dy(:,2:N-1,:) = (  d2vdy2(:,3:N,:)                    -  d2vdy2(:,1:N-2,:)) / (2*dy);
d3vdy2dy(:,    N,:) = (3*d2vdy2(:,  N,:) -4*d2vdy2(:,N-1,:) +  d2vdy2(:,  N-2,:)) / (2*dy);

duudx(1    ,:,:) = ( -uu(3  ,:,:) +4*uu( 2 ,:,:) -3*uu(1    ,:,:)) / (2*dx);
duudx(2:M-1,:,:) = (  uu(3:M,:,:)                -  uu(1:M-2,:,:)) / (2*dx);
duudx(M    ,:,:) = (3*uu(M  ,:,:) -4*uu(M-1,:,:) +  uu(  M-2,:,:)) / (2*dx);

duvdx(1    ,:,:) = ( -uv(3  ,:,:) +4*uv( 2 ,:,:) -3*uv(1    ,:,:)) / (2*dx);
duvdx(2:M-1,:,:) = (  uv(3:M,:,:)                -  uv(1:M-2,:,:)) / (2*dx);
duvdx(M    ,:,:) = (3*uv(M  ,:,:) -4*uv(M-1,:,:) +  uv(  M-2,:,:)) / (2*dx);

duvdy(:,1    ,:) = ( -uv(:,3  ,:) +4*uv(:, 2 ,:) -3*uv(:,1    ,:)) / (2*dy);
duvdy(:,2:N-1,:) = (  uv(:,3:N,:)                -  uv(:,1:N-2,:)) / (2*dy);
duvdy(:,    N,:) = (3*uv(:,  N,:) -4*uv(:,N-1,:) +  uv(:,  N-2,:)) / (2*dy);

dvvdy(:,1    ,:) = ( -vv(:,3  ,:) +4*vv(:, 2 ,:) -3*vv(:,1    ,:)) / (2*dy);
dvvdy(:,2:N-1,:) = (  vv(:,3:N,:)                -  vv(:,1:N-2,:)) / (2*dy);
dvvdy(:,    N,:) = (3*vv(:,  N,:) -4*vv(:,N-1,:) +  vv(:,  N-2,:)) / (2*dy);

d2uudx2(1    ,:,:) = (uu(3  ,:,:) -2*uu(2    ,:,:) + uu(1    ,:,:)) / (dx^2);
d2uudx2(2:M-1,:,:) = (uu(3:M,:,:) -2*uu(2:M-1,:,:) + uu(1:M-2,:,:)) / (dx^2);
d2uudx2(M    ,:,:) = (uu(  M,:,:) -2*uu(  M-1,:,:) + uu(  M-2,:,:)) / (dx^2);

d2vvdy2(:,1    ,:) = (vv(:,3  ,:) -2*vv(:,2    ,:) + vv(:,1    ,:)) / (dy^2);
d2vvdy2(:,2:N-1,:) = (vv(:,3:N,:) -2*vv(:,2:N-1,:) + vv(:,1:N-2,:)) / (dy^2);
d2vvdy2(:,    N,:) = (vv(:,  N,:) -2*vv(:,  N-1,:) + vv(:,  N-2,:)) / (dy^2);

d2udydx(1    ,:,:) = ( -dudy(3  ,:,:) +4*dudy( 2 ,:,:) -3*dudy(1    ,:,:)) / (2*dx);
d2udydx(2:M-1,:,:) = (  dudy(3:M,:,:)                  -  dudy(1:M-2,:,:)) / (2*dx);
d2udydx(M    ,:,:) = (3*dudy(M  ,:,:) -4*dudy(M-1,:,:) +  dudy(  M-2,:,:)) / (2*dx);

d2vdxdy(:,1    ,:) = ( -dvdx(:,3  ,:) +4*dvdx(:, 2 ,:) -3*dvdx(:,1    ,:)) / (2*dy);
d2vdxdy(:,2:N-1,:) = (  dvdx(:,3:N,:)                  -  dvdx(:,1:N-2,:)) / (2*dy);
d2vdxdy(:,    N,:) = (3*dvdx(:,  N,:) -4*dvdx(:,N-1,:) +  dvdx(:,  N-2,:)) / (2*dy);

%calculate wall direction, normals, and derivatives at BC pts
%still going to have problems in 1-node acute points
%will crash if narrow channels close to domain wall face outwards (ie
%must be more than 2 points existing in direction BC points towards)
for n=1:length(bclist)
    flat = 1e-6;
    
    %first point in list
%         keyboard
        xi=bclist{n}(1,1);
        yi=bclist{n}(1,2);
        s = bclist{n}(2,:) - bclist{n}(end,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(1,:) = nhat;  %store unit normal for plotting
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
        elseif nhat(1) < -flat  %east
            %use forward differences pointing left    
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
        end
        
    %middle pts
    for i=2:size(bclist{n},1)-1
        xi=bclist{n}(i,1);
        yi=bclist{n}(i,2);
        s = bclist{n}(i+1,:) - bclist{n}(i-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(i,:) = nhat;  %store unit normal for plotting
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward

        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
        end
    end %middle points
        
    %last point in list
        xi=bclist{n}(end,1);
        yi=bclist{n}(end,2);
        s = bclist{n}(1,:) - bclist{n}(end-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(end+1,:) = nhat;  %store unit normal for plotting
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
        end
end %end recalculation of derivatives on BC

%dwdz correction based on continuity
dwdz = - (dudx + dvdy);
udwdz = u.*dwdz;
vdwdz = v.*dwdz;

fprintf('%g\n',toc)
fprintf('finding extended derivatives... ')


%use new values for 1st derivatives to calculate higher time derivatives
    d2udxdt(:,:,1    ) = ( -dudx(:,:,3  ) +4*dudx(:,:, 2 ) -3*dudx(:,:,1    )) / (2*dt);
    d2udxdt(:,:,2:T-1) = (  dudx(:,:,3:T)                  -  dudx(:,:,1:T-2)) / (2*dt);
    d2udxdt(:,:,    T) = (3*dudx(:,:,  T) -4*dudx(:,:,T-1) +  dudx(:,:,  T-2)) / (2*dt);

    d2vdydt(:,:,1    ) = ( -dvdy(:,:,3  ) +4*dvdy(:,:, 2 ) -3*dvdy(:,:,1    )) / (2*dt);
    d2vdydt(:,:,2:T-1) = (  dvdy(:,:,3:T)                  -  dvdy(:,:,1:T-2)) / (2*dt);
    d2vdydt(:,:,    T) = (3*dvdy(:,:,  T) -4*dvdy(:,:,T-1) +  dvdy(:,:,  T-2)) / (2*dt);

d_udwdz_dx(1    ,:,:) = ( -udwdz(3  ,:,:) +4*udwdz( 2 ,:,:) -3*udwdz(1    ,:,:)) / (2*dx);
d_udwdz_dx(2:M-1,:,:) = (  udwdz(3:M,:,:)                   -  udwdz(1:M-2,:,:)) / (2*dx);
d_udwdz_dx(M    ,:,:) = (3*udwdz(M  ,:,:) -4*udwdz(M-1,:,:) +  udwdz(  M-2,:,:)) / (2*dx);

d_vdwdz_dy(:,1    ,:) = ( -vdwdz(:,3  ,:) +4*vdwdz(:, 2 ,:) -3*vdwdz(:,1    ,:)) / (2*dy);
d_vdwdz_dy(:,2:N-1,:) = (  vdwdz(:,3:N,:)                   -  vdwdz(:,1:N-2,:)) / (2*dy);
d_vdwdz_dy(:,    N,:) = (3*vdwdz(:,  N,:) -4*vdwdz(:,N-1,:) +  vdwdz(:,  N-2,:)) / (2*dy);

d2uvdydx(1    ,:,:) = ( -duvdy(3  ,:,:) +4*duvdy( 2 ,:,:) -3*duvdy(1    ,:,:)) / (2*dx);
d2uvdydx(2:M-1,:,:) = (  duvdy(3:M,:,:)                   -  duvdy(1:M-2,:,:)) / (2*dx);
d2uvdydx(M    ,:,:) = (3*duvdy(M  ,:,:) -4*duvdy(M-1,:,:) +  duvdy(  M-2,:,:)) / (2*dx);

d2uvdxdy(:,1    ,:) = ( -duvdx(:,3  ,:) +4*duvdx(:, 2 ,:) -3*duvdx(:,1    ,:)) / (2*dy);
d2uvdxdy(:,2:N-1,:) = (  duvdx(:,3:N,:)                   -  duvdx(:,1:N-2,:)) / (2*dy);
d2uvdxdy(:,    N,:) = (3*duvdx(:,  N,:) -4*duvdx(:,N-1,:) +  duvdx(:,  N-2,:)) / (2*dy);

%do higher order nested derivatives based on above BC derivatives
%calculate wall direction, normals, and derivatives at BC pts
%still going to have problems in 1-node acute points
%will crash if narrow channels close to domain wall face outwards (ie
%must be more than 2 points existing in direction BC points towards)
for n=1:length(bclist)
    flat = 1e-6;
    
    %first point in list
        xi=bclist{n}(1,1);
        yi=bclist{n}(1,2);
        s = bclist{n}(2,:) - bclist{n}(end,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(1,:) = nhat;  %store unit normal for plotting
        %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
        end
        
    %middle pts
    for i=2:size(bclist{n},1)-1
        xi=bclist{n}(i,1);
        yi=bclist{n}(i,2);
        s = bclist{n}(i+1,:) - bclist{n}(i-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(i,:) = nhat;  %store unit normal for plotting
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
        end
    end %middle points
        
    %last point in list
        xi=bclist{n}(end,1);
        yi=bclist{n}(end,2);
        s = bclist{n}(1,:) - bclist{n}(end-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(end+1,:) = nhat;  %store unit normal for plotting
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
        end
end %end recalculation of higher derivatives on BC

fprintf('%g\n',toc)

%%%%%%%%%%
%Calculate source terms
%  should use approach from above on BC for this, or calculate each term
%%%%%%%%%%


%Calculate N-S field function in conservative form      
if strcmp(solver,'simple')                                                  %standard BC, simplified Poisson
    f = 1/Re * (d2udx2+d2udy2) - ( dudt + u.*dudx + v.*dudy );
    g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + u.*dvdx + v.*dvdy );
    Source = -( dudx.*dudx + dudy.*dvdx ...
               +dvdx.*dudy + dvdy.*dvdy );
   
elseif strcmp(solver,'simple2')                                             %standard BC, simplified Poisson+dwdz
    f = 1/Re * (d2udx2+d2udy2) - ( dudt + u.*dudx + v.*dudy );
    g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + u.*dvdx + v.*dvdy );
    Source = -( dudx.*dudx + dudy.*dvdx ...
               +dvdx.*dudy + dvdy.*dvdy ...
               +dwdz.*dwdz );   
elseif strcmp(solver,'standard')                                            %standard BC, full standard Poisson
    f = 1/Re * (d2udx2+d2udy2) - ( dudt + u.*dudx + v.*dudy );
    g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + u.*dvdx + v.*dvdy );
    Source = 1/Re*( d3udx2dx + d3udy2dx + d3vdx2dy + d3vdy2dy) ...
                - ( d2udxdt + d2vdydt ...
                   +dudx.*dudx + dudy.*dvdx ...
                   +dvdx.*dudy + dvdy.*dvdy ...
                   +u.*d2udx2 + v.*d2udydx + u.*d2vdxdy + v.*d2vdy2 );
elseif strcmp(solver,'standard2')                                           %conservative convc. BC, full standard Poisson
    f = 1/Re * (d2udx2+d2udy2) - ( dudt + duudx + duvdy + udwdz);
    g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + duvdx + dvvdy + vdwdz);
    Source = 1/Re*( d3udx2dx + d3udy2dx + d3vdx2dy + d3vdy2dy) ...
                - ( d2udxdt + d2vdydt ...
                   +dudx.*dudx + dudy.*dvdx ...
                   +dvdx.*dudy + dvdy.*dvdy ...
                   +u.*d2udx2 + v.*d2udydx + u.*d2vdxdy + v.*d2vdy2 );
elseif strcmp(solver,'conservative')                                        %conservative convection BC and Poisson
    f = 1/Re * (d2udx2+d2udy2) - ( dudt + duudx + duvdy + udwdz);
    g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + duvdx + dvvdy + vdwdz);
    Source = 1/Re*( d3udx2dx + d3udy2dx + d3vdx2dy + d3vdy2dy) ...
                - ( d2udxdt + d2vdydt ...
                   +d2uudx2  + d2uvdydx + d_udwdz_dx ...
                   +d2uvdxdy + d2vvdy2  + d_vdwdz_dy );
else
    error('unknown solver');
end

          
%keyboard

% clear dudt dvdt dudx dvdx dudy dvdy d2udx2 d2udy2 d2vdx2 d2vdy2
% clear dfdx dgdy uu uv vv udwdz vdwdz
% clear d2udydx d3udx2dx d3udy2dx d2vdxdy d3vdx2dy d3vdy2dy d2udxdt d2vdydt
% clear d2uudx2 d2uvdydx d_udwdz_dx d2uvdxdy d2vvdy2 d_vdwdz_dy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use SOC for middle of region
% the problem has been discritized so that

%default discretization
Ap(1:NI,1:NJ)   = -2*(1/dx^2+1/dy^2);
Ae(1:NI,1:NJ)   = 1/dx^2;
Aw(1:NI,1:NJ)   = 1/dx^2;
An(1:NI,1:NJ)   = 1/dy^2;
As(1:NI,1:NJ)   = 1/dy^2;
 S(1:NI,1:NJ,:) = Source(1:NI,1:NJ,:);


%calculate wall direction, normals, and BC
for n=1:length(bclist)
    
    %first point in list
        xi=bclist{n}(1,1);
        yi=bclist{n}(1,2);
        s = bclist{n}(2,:) - bclist{n}(end,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        
        if nhat(1) > 0          %west
            if nhat(2) > 0      % south
                Ap(xi, yi)   = -1/dx*nhat(1) - 1/dy*nhat(2);
                Ae(xi, yi)   =  1/dx*nhat(1);
                Aw(xi, yi)   =     0;
                An(xi, yi)   =                 1/dy*nhat(2);
                As(xi, yi)   =     0;
            else                % north
                Ap(xi, yi)   = -1/dx*nhat(1) + 1/dy*nhat(2);
                Ae(xi, yi)   =  1/dx*nhat(1);
                Aw(xi, yi)   =     0;
                An(xi, yi)   =     0;
                As(xi, yi)   =               - 1/dy*nhat(2);   
            end
            
        else                    %east
            if nhat(2) > 0      % south
                Ap(xi, yi)   =  1/dx*nhat(1) - 1/dy*nhat(2);
                Ae(xi, yi)   =     0;
                Aw(xi, yi)   = -1/dx*nhat(1);
                An(xi, yi)   =                 1/dy*nhat(2);
                As(xi, yi)   =     0;
            else                % north
                Ap(xi, yi)   =  1/dx*nhat(1) + 1/dy*nhat(2);
                Ae(xi, yi)   =     0;
                Aw(xi, yi)   = -1/dx*nhat(1);
                An(xi, yi)   =     0;
                As(xi, yi)   =                -1/dy*nhat(2);   
            end
        end    
        S(xi, yi,:)  =  f(xi, yi,:)*nhat(1) + g(xi,yi,:)*nhat(2);
        
    %middle pts
    for i=2:size(bclist{n},1)-1
        xi=bclist{n}(i,1);
        yi=bclist{n}(i,2);
        s = bclist{n}(i+1,:) - bclist{n}(i-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        
        if nhat(1) > 0          %west
            if nhat(2) > 0      % south
                Ap(xi, yi)   = -1/dx*nhat(1) - 1/dy*nhat(2);
                Ae(xi, yi)   =  1/dx*nhat(1);
                Aw(xi, yi)   =     0;
                An(xi, yi)   =                 1/dy*nhat(2);
                As(xi, yi)   =     0;

            else                % north
                Ap(xi, yi)   = -1/dx*nhat(1) + 1/dy*nhat(2);
                Ae(xi, yi)   =  1/dx*nhat(1);
                Aw(xi, yi)   =     0;
                An(xi, yi)   =     0;
                As(xi, yi)   =               - 1/dy*nhat(2);   
            end
            
        else                    %east
            if nhat(2) > 0      % south
                Ap(xi, yi)   =  1/dx*nhat(1) - 1/dy*nhat(2);
                Ae(xi, yi)   =     0;
                Aw(xi, yi)   = -1/dx*nhat(1);
                An(xi, yi)   =                 1/dy*nhat(2);
                As(xi, yi)   =     0;
            else                % north
                Ap(xi, yi)   =  1/dx*nhat(1) + 1/dy*nhat(2);
                Ae(xi, yi)   =     0;
                Aw(xi, yi)   = -1/dx*nhat(1);
                An(xi, yi)   =     0;
                As(xi, yi)   =                -1/dy*nhat(2);   
            end
        end
        S(xi, yi,:)  =  f(xi, yi,:)*nhat(1) + g(xi,yi,:)*nhat(2);       
    end

    %last point in list
        xi=bclist{n}(end,1);
        yi=bclist{n}(end,2);
        s = bclist{n}(1,:) - bclist{n}(end-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        
        if nhat(1) > 0          %west
            if nhat(2) > 0      % south
                Ap(xi, yi)   = -1/dx*nhat(1) - 1/dy*nhat(2);
                Ae(xi, yi)   =  1/dx*nhat(1);
                Aw(xi, yi)   =     0;
                An(xi, yi)   =                 1/dy*nhat(2);
                As(xi, yi)   =     0;

            else                % north
                Ap(xi, yi)   = -1/dx*nhat(1) + 1/dy*nhat(2);
                Ae(xi, yi)   =  1/dx*nhat(1);
                Aw(xi, yi)   =     0;
                An(xi, yi)   =     0;
                As(xi, yi)   =               - 1/dy*nhat(2);   
            end
            
        else                    %east
            if nhat(2) > 0      % south
                Ap(xi, yi)   =  1/dx*nhat(1) - 1/dy*nhat(2);
                Ae(xi, yi)   =     0;
                Aw(xi, yi)   = -1/dx*nhat(1);
                An(xi, yi)   =                 1/dy*nhat(2);
                As(xi, yi)   =     0;
            else                % north
                Ap(xi, yi)   =  1/dx*nhat(1) + 1/dy*nhat(2);
                Ae(xi, yi)   =     0;
                Aw(xi, yi)   = -1/dx*nhat(1);
                An(xi, yi)   =     0;
                As(xi, yi)   =                -1/dy*nhat(2);   
            end
        end
        S(xi, yi,:)  =  f(xi, yi,:)*nhat(1) + g(xi,yi,:)*nhat(2);
end %for n=1:length(bclist) 

% %Bad idea - convergence is terrible - just cheat and shift the P at Px,Py
% turn it back on after conversation with Osama
if nargin>=9
    %Fixed pressure measurement at a point
    Ap(Px,Py)   = 1;
    Ae(Px,Py)   = 0;
    Aw(Px,Py)   = 0;
    An(Px,Py)   = 0;
    As(Px,Py)   = 0;
     S(Px,Py,:) = Pref;
end

%convert coefficient lists to a sparse matrix
AA = spdiags( [reshape(As,NI*NJ,1), ...
               reshape(Aw,NI*NJ,1), ...
               reshape(Ap,NI*NJ,1), ...
               reshape(Ae,NI*NJ,1), ...
               reshape(An,NI*NJ,1)] , [+NI,+1,0,-1,-NI], NI*NJ, NJ*NI ).';

ktot = 0;

pause(0.001)

%solve for pressures
fprintf('     t       k       deltaP  flag        time\n')
fprintf('  ----  ------  -----------  ----  ----------\n')

for t=1:NT
    
    [phi]=AA\reshape(S(:,:,t),NI*NJ,1);  flag=0;relres=0;iter=[1,1]; resvec=0;

    %store final solution into pressure table
    p(:,:,t) = reshape(phi,NI,NJ);
%     keyboard
    %use an offset to peg the reference pressure over time... 
    %is this okay to do?
    if nargin>=9
        deltaP = Pref(t) - p(Px,Py,t);
        p(:,:,t) = p(:,:,t) + deltaP * ones(NI,NJ);
    else
        deltaP = 0;
    end

    %convergence information
    %fprintf('  %4i  %6i  %11.4g  %4i  %10.3f\n',t,(iter(1)-1)*restarts+iter(2),relres,flag,toc)
    fprintf('  %4i  %6i  %11.4g  %4i  %10.3f\n',t,(iter(1)-1)*restarts+iter(2),deltaP,flag,toc)
    %plot current result for testing
%     figure(944),mesh(x,y,p(:,:,t)),pause(0.0001),colorbar
%     figure(945),imagesc([x(1,1),x(end,end)],[y(1,1),y(end,end)],p(:,:,t)),pause(0.0001),axis xy, colorbar
    %figure(945),contourf(x,y,p(:,:,t),20),pause(0.0001)

%     %plot convergence
%     normS = norm(reshape(S(:,:,t),NI*NJ,1));
%     figure(942),hold on,semilogy(ktot+[0:length(resvec)-1],resvec/normS),title(['residual for timestep t=',num2str(t)])
%     figure(943),semilogy(resvec/normS),title(['residual for timestep t=',num2str(t)])
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
