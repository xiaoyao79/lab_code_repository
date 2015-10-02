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
p = zeros(NI,NJ,NT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the source term 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_size = size(u);   %u and v should be the same size
M = u_size(1);      %number of points in x
N = u_size(2);      %number of points in y
T = u_size(3);      %number of points in t

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
dvdx = compactdiff(v,dx,1);
dudy = compactdiff(u,dy,2);
dvdy = compactdiff(v,dy,2);

d2udx2 = compactdiff2(u,dx,1);
d2udy2 = compactdiff2(u,dy,2);
d2vdx2 = compactdiff2(v,dx,1);
d2vdy2 = compactdiff2(v,dy,2);

duudx = compactdiff(uu,dx,1);
duvdx = compactdiff(uv,dx,1);
duvdy = compactdiff(uv,dy,2);
dvvdy = compactdiff(vv,dy,2);

d2udydx = compactdiff(dudy,dx,1);
d2vdxdy = compactdiff(dvdx,dy,2);

%dwdz correction based on continuity
dwdz = - (dudx + dvdy);
udwdz = u.*dwdz;
vdwdz = v.*dwdz;


d_dwdz_dx = compactdiff(dwdz,dx,1);
d_dwdz_dy = compactdiff(dwdz,dy,2);

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

%can probably prune the list of derivatives for line integrals - many terms
%are for the source term only
clear dudt dvdt dudx dvdx dudy dvdy d2udx2 d2udy2 d2vdx2 d2vdy2
clear dfdx dgdy uu uv vv udwdz vdwdz
clear d2udydx d3udx2dx d3udy2dx d2vdxdy d3vdx2dy d3vdy2dy d2udxdt d2vdydt
clear d2uudx2 d2uvdydx d_udwdz_dx d2uvdxdy d2vvdy2 d_vdwdz_dy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ktot = 0;
pause(0.001)
%solve for pressures
fprintf('     t       k        resid  flag        time\n')
fprintf('  ----  ------  -----------  ----  ----------\n')

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
    
    flag=0;relres=0;iter=[1,1];,resvec=0;

    %store final solution into pressure table
    p(:,:,t) = (p1+p2+p3+p4+p5+p6+p7+p8)/8;
    
    if nargin>=9
        %pause
        deltaP = Pref(t) - p(Px,Py,t);
        p(:,:,t) = p(:,:,t) + deltaP * ones(NI,NJ);
    end

    %convergence information
    fprintf('  %4i  %6i  %11.4g  %4i  %10.3f\n',t,(iter(1)-1)*restarts+iter(2),relres,flag,toc)
    pause(1e-4)
       
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
