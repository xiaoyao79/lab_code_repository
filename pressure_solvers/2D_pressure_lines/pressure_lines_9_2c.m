function [p]=pressure_lines_9_2c(x,y,u,v,Re,dt,solver,bcmask,Pref,Px,Py)
%function [P]=pressure_lines_9_2(X,Y,U,V,Re,DT,[SOLVER],[BCMASK],[PREF, PX, PY])
%
% Returns P(x,y,t) from U(x,y,t) and V(x,y,t) using a multi-directional 
% Guass-Seidel style of line integration on pressure terms of the Navier-
% Stokes Equations.  DT is the non-dimensional time step between velocity 
% snapshots.  By default, a conservative formulation.
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
% not have to be the same, and are determined from X and Y.  Note that the
% positive X direction must increment along the 1st index, and the positive
% Y direction must increment along the second index for 
%
% All variables should be normalized by:
%   U = U'/U0, V = V'/U0, X = X'/L, Y = Y'/L
%   t = T'/T , T=L/U0, Re = L*U0/nu, DT=dt'/T
%   P = P'/P0, P0 = rho*U0^2
%
% The choice of normalization must be consistent throughout the inputs. 
% For example, choosing U0=1 and L=1 (physical units on length and speed) 
% results in T=1 (and therefore DT=dt'), but Re = 1/nu.
%
% BCMASK should be a binary mask indicating which vector locations are to
% be used in the calculation of pressure.  It is ordered to match U and V.
% Degenerate shapes (single isolated vector locations), single vector 
% lines, and complicated corners near the edge of the X-Y domain may cause 
% problems.  
% BCMASK(x,y) results in the same boundary conditions being used at each
% time step.
% BCMASK(x,y,t) generates a new boundary condition for each time step, and
% must be the same size as U and V.
%
% If BCMASK is not defined, the default is to use the domain walls as BC's
%
% PREF is a normalized pressure measured at index pt [PX,PY]
% If PREF is a single value, it is assumed that the pressure at [PX,PY] 
% is constant for all time steps.  Otherwise it must have the same number 
% of time steps as U and V. PREF is optional, but if not defined, the
% pressure will float between time steps and cannot be compared.

% obsolete syntax:
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
% v7: converting to use omni-directional integration on larger outer domain
% as per Liu and Katz
% v8: adding iteration to boundary condition loop pre-line integration
% v8_1: change to radius of 1.5 (from 1.2), and NR=2*(NI+NJ)
%       added performance enhancements, inlined linspace() and s2i(), about
%       30% speed increase
% v9_1: switch to assuming  periodic in time - ie U(t=T+1)=U(t=1);
% v9_2: allow time-varying BC's (duplicated from v8_2)
% v9_2a: keep track of convergence rate as well as residual, plot both
%        add termination criterion for convergence rate
% v9_2b: also check calculated pressure gradient as compared to source term

warning('off','MATLAB:Axes:NegativeDataInLogAxis');

% try
DIFFSCHEME = 'soc';
RSIZE = 1.5;  %factor for increasing radius of integrating circle
%NR = 200;

%figure(1234),quiver(x,y,u(:,:,1),v(:,:,1)),pause

NI = size(x,1);
NJ = size(y,2);
NT = size(u,3);

%NR      = 2*(NI+NJ);
NR      = 2*(NI+NJ)/2;   %number of points on perimeter of virtual boundary
QMAX    = 50;            %maximum number of iterations per snapshot (overrides convergence)
tol     = 5e-3;          %convergence criterion on residuals
tolcon  = 1-1e-2;          %convergence criterion on convergence rate

%dynamic boundary conditions based on function arguments
%Neumann-type pressure gradient BC's
if nargin>=7    %required args
    if nargin>=8    %includes bcmask
        if nargin>=9   %includes at least 1 pressure arg
            if nargin<11    %but doesn't include all
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
%         %disp('BCLIST is defined')
%         if iscell(bclist)
%             %disp('BCLIST is a cell')
%             if size(bclist) == 0;
%                 disp('BCLIST must have at least 1 element')  
%                 error('BCLIST must have the following form: { [x1,y1;...],[x1,y2;...], ... }')
%             end 
%         else %can't use this BCLIST
%             disp('BCLIST is not a cell')                        
%             error('BCLIST must have the following form: { [x1,y1;...],[x1,y2;...], ... }')
%         end

    else %didn't include any bc's
        disp('BCLIST is not defined, using default boundary conditions')
            %define BCLIST to use the domain edges as BC's
%             bclist{1} = [[      (1:NI-1)' , NJ*ones(NI-1,1) ];...    %north
%                          [NI*ones(NJ-1,1) ,       (NJ:-1:2)'];...    %east  
%                          [      (NI:-1:2)',    ones(NI-1,1) ];...    %south
%                          [   ones(NJ-1,1),         (1:NJ-1)']];     %west
            %define BCMASK to use the domain edges as BC's
            bcmask = ones(NI,NJ);         
    end
%     keyboard
    %check dimensions of BCMASK and replicate BCMASK, if needed
    if ndims(bcmask)==2 & size(bcmask)==size(x)
        bcmask = repmat(logical(bcmask),[1,1,NT]);  %truncate bcmask to logical
    elseif ndims(bcmask)==3 & size(bcmask)==size(u)
        bcmask = logical(bcmask);                   %truncate bcmask to logical
    else
        error('If size(U)==[NX,NY,NY], then size(BCMASK) must be either [NX,NY,1] or [NX,NY,NT]')
    end
    
    %fill bclist from bcmask, reorder any inner boundaries
    bclist = cell(1,NT);
    for k=1:NT
        [bclist{k},bcL,bcN,bcA] = bwboundaries(bcmask(:,:,k),8,'holes');  %bcL is labels, bcN is #/objects, bcA is adjacency matrix
        
        for n=1:bcN %first bcN entries are correctly ordered, truncate duplicate endpoint
            %bclist{k}{n} = flipdim(flipdim(bclist{k}{n}(1:end-1,:),1),2);
            bclist{k}{n} = bclist{k}{n}(1:end-1,:);
        end
        for n=bcN+1:length(bclist{k}) %remaining entries are reverse ordered (holes), truncate duplicate endpoint
            %bclist{k}{n} = flipdim(bclist{k}{n}(1:end-1,:),2);
            bclist{k}{n} = flipdim(bclist{k}{n}(1:end-1,:),1);
        end
    end
    
else
    %flowtype = 'viscous'
    solver = 'conservative';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
tic

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

% % storage for source term
% S = zeros(NI,NJ,NT);
% % phi(i,j,k) where i is x-index, j is y-index, and k=1,2,3 for time-stepping
% %discretization defined on the following neighboring points
% AA = sparse(NI*NJ,NI*NJ);
% An = zeros(NI,NJ);  
% As = zeros(NI,NJ);
% Aw = zeros(NI,NJ);
% Ae = zeros(NI,NJ);
% Ap = zeros(NI,NJ);
% %store pressure in phi during solution
% phi = zeros(NI*NJ,1);
%source term plus LHS terms not used in tridiagonal solution
% RHSi = zeros(NI,1);
% RHSj = zeros(NJ,1);
%storage for final pressure solutions
p = zeros(NI,NJ,NT);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the source term 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_size = size(u);   %u and v should be the same size
% M = u_size(1);      %number of points in x
% N = u_size(2);      %number of points in y
T = u_size(3);      %number of points in t - same as NT above

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

%change time derivatives to assume periodic
if NT<7
    %dudt(:,:,1    ) = ( -u(:,:,3  ) +4*u(:,:, 2 ) -3*u(:,:,1    )) / (2*dt);
    dudt(:,:,1    ) = (  u(:,:,2  )               -  u(:,:,T    )) / (2*dt);
    dudt(:,:,2:T-1) = (  u(:,:,3:T)               -  u(:,:,1:T-2)) / (2*dt);
    dudt(:,:,    T) = (  u(:,:,  1)               -  u(:,:,  T-1)) / (2*dt);
    %dudt(:,:,    T) = (3*u(:,:,  T) -4*u(:,:,T-1) +  u(:,:,  T-2)) / (2*dt);
    
    %dvdt(:,:,1    ) = ( -v(:,:,3  ) +4*v(:,:, 2 ) -3*v(:,:,1    )) / (2*dt);
    dvdt(:,:,1    ) = (  v(:,:,2  )               -  v(:,:,T    )) / (2*dt);
    dvdt(:,:,2:T-1) = (  v(:,:,3:T)               -  v(:,:,1:T-2)) / (2*dt);
    dvdt(:,:,    T) = (  v(:,:,  1)               -  v(:,:,  T-1)) / (2*dt);
    %dvdt(:,:,    T) = (3*v(:,:,  T) -4*v(:,:,T-1) +  v(:,:,  T-2)) / (2*dt);
else
    if strcmp(DIFFSCHEME,'compact')
        dudt = compactdiff(u,dt,3);
        dvdt = compactdiff(v,dt,3);
    elseif strcmp(DIFFSCHEME,'soc')
        %dudt = socdiff(u,dt,3);
        %dvdt = socdiff(v,dt,3);
%         dudt(:,:,1    ) = (  u(:,:,2  )               -  u(:,:,T    )) / (2*dt);
%         dudt(:,:,2:T-1) = (  u(:,:,3:T)               -  u(:,:,1:T-2)) / (2*dt);
%         dudt(:,:,    T) = (  u(:,:,  1)               -  u(:,:,  T-1)) / (2*dt);
%         dvdt(:,:,1    ) = (  v(:,:,2  )               -  v(:,:,T    )) / (2*dt);
%         dvdt(:,:,2:T-1) = (  v(:,:,3:T)               -  v(:,:,1:T-2)) / (2*dt);
%         dvdt(:,:,    T) = (  v(:,:,  1)               -  v(:,:,  T-1)) / (2*dt);

        %need to loop first and last 2 frames in order to allow proper
        %handling of periodicity - need up to 2 frames forward or back
        

        temp = socdiff_bc( cat(3,      u(:,:,end-1),      u(:,:,end),      u,      u(:,:,1),      u(:,:,2)), dt, 3, ...
                           cat(3, bcmask(:,:,end-1), bcmask(:,:,end), bcmask, bcmask(:,:,1), bcmask(:,:,2)) );
        dudt = temp(:,:,3:end-2);
        temp = socdiff_bc( cat(3,      v(:,:,end-1),      v(:,:,end),      v,      v(:,:,1),      v(:,:,3)), dt, 3, ...
                           cat(3, bcmask(:,:,end-1), bcmask(:,:,end), bcmask, bcmask(:,:,1), bcmask(:,:,2)) );
        dvdt = temp(:,:,3:end-2);

    end
end

if strcmp(DIFFSCHEME,'compact')
    dudx   = compactdiff(u,dx,1);
    dvdx   = compactdiff(v,dx,1);
    dudy   = compactdiff(u,dy,2);
    dvdy   = compactdiff(v,dy,2);
    d2udx2 = compactdiff2(u,dx,1);
    d2udy2 = compactdiff2(u,dy,2);
    d2vdx2 = compactdiff2(v,dx,1);
    d2vdy2 = compactdiff2(v,dy,2);

    % d3udx2dx = compactdiff(d2udx2,dx,1);%
    % d3udy2dx = compactdiff(d2udy2,dx,1);
    % d3vdx2dy = compactdiff(d2vdx2,dy,2);
    % d3vdy2dy = compactdiff(d2vdy2,dy,2);

    duudx = compactdiff(uu,dx,1);
    duvdx = compactdiff(uv,dx,1);
    duvdy = compactdiff(uv,dy,2);
    dvvdy = compactdiff(vv,dy,2);
    % d2uudx2 = compactdiff2(uu,dx,1);
    % d2vvdy2 = compactdiff2(vv,dy,2);

    d2udydx = compactdiff(dudy,dx,1);
    d2vdxdy = compactdiff(dvdx,dy,2);
elseif strcmp(DIFFSCHEME,'soc')
    dudx   = socdiff(u,dx,1);
    dvdx   = socdiff(v,dx,1);
    dudy   = socdiff(u,dy,2);
    dvdy   = socdiff(v,dy,2);
    d2udx2 = socdiff2(u,dx,1);
    d2udy2 = socdiff2(u,dy,2);
    d2vdx2 = socdiff2(v,dx,1);
    d2vdy2 = socdiff2(v,dy,2);

    % d3udx2dx = socdiff(d2udx2,dx,1);%
    % d3udy2dx = socdiff(d2udy2,dx,1);
    % d3vdx2dy = socdiff(d2vdx2,dy,2);
    % d3vdy2dy = socdiff(d2vdy2,dy,2);

    duudx = socdiff(uu,dx,1);
    duvdx = socdiff(uv,dx,1);
    duvdy = socdiff(uv,dy,2);
    dvvdy = socdiff(vv,dy,2);
    % d2uudx2 = socdiff2(uu,dx,1);
    % d2vvdy2 = socdiff2(vv,dy,2);

    d2udydx = socdiff(dudy,dx,1);
    d2vdxdy = socdiff(dvdx,dy,2);
end

%% correct BC for lower order derivatives
%***THESE SECTIONS MIGHT BE BETTER HANDLED SWITCHING TO socdiff_bc ABOVE***
%calculate wall direction, normals, and derivatives at BC pts
%still going to have problems in 1-node acute points
%will crash if narrow channels close to domain wall face outwards (ie
%must be more than 2 points existing in direction BC points towards)
for t=1:NT
for n=1:1
    flat = 1e-6;
    %disp(' xi  yi       s1       s2    nhat1    nhat2  N/S  E/W')
    
    %first point in list
        xi=bclist{t}{n}(1,1);
        yi=bclist{t}{n}(1,2);
%         keyboard
        if size(bclist{t}{n},1)==1
            s = [1 1];
        else
            s = bclist{t}{n}(2,:) - bclist{t}{n}(end,:);
        end
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(1,:) = nhat;  %store unit normal for plotting
        %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward
            %fprintf('     N')
            dudy  (xi, yi,:)    = (3*u(xi, yi,:) -4*u(xi, yi-1,:) + u(xi, yi-2,:)) / (2*dy);
            dvdy  (xi, yi,:)    = (3*v(xi, yi,:) -4*v(xi, yi-1,:) + v(xi, yi-2,:)) / (2*dy);
            d2udy2(xi, yi,:)    = (  u(xi, yi,:) -2*u(xi, yi-1,:) + u(xi, yi-2,:)) / (dy^2);
            d2vdy2(xi, yi,:)    = (  v(xi, yi,:) -2*v(xi, yi-1,:) + v(xi, yi-2,:)) / (dy^2);
            duvdy  (xi, yi,:)    = (3*uv(xi, yi,:) -4*uv(xi, yi-1,:) + uv(xi, yi-2,:)) / (2*dy);
            dvvdy  (xi, yi,:)    = (3*vv(xi, yi,:) -4*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (2*dy);
            %d2vvdy2(xi, yi,:)    = (  vv(xi, yi,:) -2*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (dy^2);
            
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
            %fprintf('     S')
            dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dy);
            dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dy);
            d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
            d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
            duvdy  (xi, yi,:) = (-uv(xi, yi+2,:) +4*uv(xi, yi+1,:) -3*uv(xi, yi,:)) / (2*dy);
            dvvdy  (xi, yi,:) = (-vv(xi, yi+2,:) +4*vv(xi, yi+1,:) -3*vv(xi, yi,:)) / (2*dy);
            %d2vvdy2(xi, yi,:) = ( vv(xi, yi+2,:) -2*vv(xi, yi+1,:) +  vv(xi, yi,:)) / (dy^2);

        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
            %fprintf('     0')
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
            %fprintf('     W\n')
            dudx  (xi ,yi,:) = (-u(xi+2, yi,:) +4*u(xi+1, yi,:) -3*u(xi, yi,:)) / (2*dx);
            dvdx  (xi ,yi,:) = (-v(xi+2, yi,:) +4*v(xi+1, yi,:) -3*v(xi, yi,:)) / (2*dx);
            d2udx2(xi ,yi,:) = ( u(xi+2, yi,:) -2*u(xi+1, yi,:) +  u(xi, yi,:)) / (dx^2);
            d2vdx2(xi ,yi,:) = ( v(xi+2, yi,:) -2*v(xi+1, yi,:) +  v(xi, yi,:)) / (dx^2);
            duudx  (xi ,yi,:) = (-uu(xi+2, yi,:) +4*uu(xi+1, yi,:) -3*uu(xi, yi,:)) / (2*dx);
            duvdx  (xi ,yi,:) = (-uv(xi+2, yi,:) +4*uv(xi+1, yi,:) -3*uv(xi, yi,:)) / (2*dx);
            %d2uudx2(xi ,yi,:) = ( uu(xi+2, yi,:) -2*uu(xi+1, yi,:) +  uu(xi, yi,:)) / (dx^2);

        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
            %fprintf('     E\n')
            dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
            dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
            d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
            d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
            duudx  (xi, yi,:) = (3*uu(xi, yi,:) -4*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (2*dx);
            duvdx  (xi, yi,:) = (3*uv(xi, yi,:) -4*uv(xi-1, yi,:) + uv(xi-2, yi,:)) / (2*dx);
            %d2uudx2(xi, yi,:) = (  uu(xi, yi,:) -2*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (dx^2);
            
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
            %fprintf('     0\n')
        end
        %pause
        
    %middle pts
    for i=2:size(bclist{t}{n},1)-1
        xi=bclist{t}{n}(i,1);
        yi=bclist{t}{n}(i,2);
        s = bclist{t}{n}(i+1,:) - bclist{t}{n}(i-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(i,:) = nhat;  %store unit normal for plotting
        %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward
            %fprintf('     N')
            dudy  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi, yi-1,:) + u(xi, yi-2,:)) / (2*dy);
            dvdy  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi, yi-1,:) + v(xi, yi-2,:)) / (2*dy);
            d2udy2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi, yi-1,:) + u(xi, yi-2,:)) / (dy^2);
            d2vdy2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi, yi-1,:) + v(xi, yi-2,:)) / (dy^2);
            duvdy  (xi, yi,:)    = (3*uv(xi, yi,:) -4*uv(xi, yi-1,:) + uv(xi, yi-2,:)) / (2*dy);
            dvvdy  (xi, yi,:)    = (3*vv(xi, yi,:) -4*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (2*dy);
            %d2vvdy2(xi, yi,:)    = (  vv(xi, yi,:) -2*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (dy^2);
            
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
            %fprintf('     S')
            dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dy);
            dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dy);
            d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
            d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
            duvdy  (xi, yi,:) = (-uv(xi, yi+2,:) +4*uv(xi, yi+1,:) -3*uv(xi, yi,:)) / (2*dy);
            dvvdy  (xi, yi,:) = (-vv(xi, yi+2,:) +4*vv(xi, yi+1,:) -3*vv(xi, yi,:)) / (2*dy);
            %d2vvdy2(xi, yi,:) = ( vv(xi, yi+2,:) -2*vv(xi, yi+1,:) +  vv(xi, yi,:)) / (dy^2);
        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
            %fprintf('     0')
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
            %fprintf('     W\n')
            dudx  (xi ,yi,:) = (-u(xi+2, yi,:) +4*u(xi+1, yi,:) -3*u(xi, yi,:)) / (2*dx);
            dvdx  (xi ,yi,:) = (-v(xi+2, yi,:) +4*v(xi+1, yi,:) -3*v(xi, yi,:)) / (2*dx);
            d2udx2(xi ,yi,:) = ( u(xi+2, yi,:) -2*u(xi+1, yi,:) +  u(xi, yi,:)) / (dx^2);
            d2vdx2(xi ,yi,:) = ( v(xi+2, yi,:) -2*v(xi+1, yi,:) +  v(xi, yi,:)) / (dx^2);
            duudx  (xi ,yi,:) = (-uu(xi+2, yi,:) +4*uu(xi+1, yi,:) -3*uu(xi, yi,:)) / (2*dx);
            duvdx  (xi ,yi,:) = (-uv(xi+2, yi,:) +4*uv(xi+1, yi,:) -3*uv(xi, yi,:)) / (2*dx);
            %d2uudx2(xi ,yi,:) = ( uu(xi+2, yi,:) -2*uu(xi+1, yi,:) +  uu(xi, yi,:)) / (dx^2);

        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
            %fprintf('     E\n')
            dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
            dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
            d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
            d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
            duudx  (xi, yi,:) = (3*uu(xi, yi,:) -4*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (2*dx);
            duvdx  (xi, yi,:) = (3*uv(xi, yi,:) -4*uv(xi-1, yi,:) + uv(xi-2, yi,:)) / (2*dx);
            %d2uudx2(xi, yi,:) = (  uu(xi, yi,:) -2*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (dx^2);
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
            %fprintf('     0\n')
        end
        %pause
    end %middle points
        
    %last point in list
        xi=bclist{t}{n}(end,1);
        yi=bclist{t}{n}(end,2);
        s = bclist{t}{n}(1,:) - bclist{t}{n}(end-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(end+1,:) = nhat;  %store unit normal for plotting
        %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)
        
        if nhat(2) < -flat      % north
            %use forward differences pointing downward
            %fprintf('     N')
            dudy  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi, yi-1,:) + u(xi, yi-2,:)) / (2*dy);
            dvdy  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi, yi-1,:) + v(xi, yi-2,:)) / (2*dy);
            d2udy2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi, yi-1,:) + u(xi, yi-2,:)) / (dy^2);
            d2vdy2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi, yi-1,:) + v(xi, yi-2,:)) / (dy^2);
            duvdy  (xi, yi,:)    = (3*uv(xi, yi,:) -4*uv(xi, yi-1,:) + uv(xi, yi-2,:)) / (2*dy);
            dvvdy  (xi, yi,:)    = (3*vv(xi, yi,:) -4*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (2*dy);
            %d2vvdy2(xi, yi,:)    = (  vv(xi, yi,:) -2*vv(xi, yi-1,:) + vv(xi, yi-2,:)) / (dy^2);
           
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
            %fprintf('     S')
            dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dy);
            dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dy);
            d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
            d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
            duvdy  (xi, yi,:) = (-uv(xi, yi+2,:) +4*uv(xi, yi+1,:) -3*uv(xi, yi,:)) / (2*dy);
            dvvdy  (xi, yi,:) = (-vv(xi, yi+2,:) +4*vv(xi, yi+1,:) -3*vv(xi, yi,:)) / (2*dy);
            %d2vvdy2(xi, yi,:) = ( vv(xi, yi+2,:) -2*vv(xi, yi+1,:) +  vv(xi, yi,:)) / (dy^2);
        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
            %fprintf('     0')
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
            %fprintf('     W\n')
            dudx  (xi ,yi,:) = (-u(xi+2, yi,:) +4*u(xi+1, yi,:) -3*u(xi, yi,:)) / (2*dx);
            dvdx  (xi ,yi,:) = (-v(xi+2, yi,:) +4*v(xi+1, yi,:) -3*v(xi, yi,:)) / (2*dx);
            d2udx2(xi ,yi,:) = ( u(xi+2, yi,:) -2*u(xi+1, yi,:) +  u(xi, yi,:)) / (dx^2);
            d2vdx2(xi ,yi,:) = ( v(xi+2, yi,:) -2*v(xi+1, yi,:) +  v(xi, yi,:)) / (dx^2);
            duudx  (xi ,yi,:) = (-uu(xi+2, yi,:) +4*uu(xi+1, yi,:) -3*uu(xi, yi,:)) / (2*dx);
            duvdx  (xi ,yi,:) = (-uv(xi+2, yi,:) +4*uv(xi+1, yi,:) -3*uv(xi, yi,:)) / (2*dx);
            %d2uudx2(xi ,yi,:) = ( uu(xi+2, yi,:) -2*uu(xi+1, yi,:) +  uu(xi, yi,:)) / (dx^2);

        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
            %fprintf('     E\n')
            dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
            dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
            d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
            d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
            duudx  (xi, yi,:) = (3*uu(xi, yi,:) -4*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (2*dx);
            duvdx  (xi, yi,:) = (3*uv(xi, yi,:) -4*uv(xi-1, yi,:) + uv(xi-2, yi,:)) / (2*dx);
            %d2uudx2(xi, yi,:) = (  uu(xi, yi,:) -2*uu(xi-1, yi,:) + uu(xi-2, yi,:)) / (dx^2);
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
            %fprintf('     0\n')
        end
        %pause
%     figure(2345)
%     hold on
%     plot(bclist{t}{n}(:,1),bclist{t}{n}(:,2))
%     quiver(bclist{t}{n}(:,1),bclist{t}{n}(:,2),nmlist{n}(:,1),nmlist{n}(:,2))
%     title(num2str(n))
%     pause
%     keyboard
end %end recalculation of derivatives on BC
% figure(2345)
% hold off
end

%% higher order derivatives
%dwdz correction based on continuity
dwdz = - (dudx + dvdy);
udwdz = u.*dwdz;
vdwdz = v.*dwdz;

if strcmp(solver,'conservative2')
    if strcmp(DIFFSCHEME,'compact')

        % %use new values for 1st derivatives to calculate higher time derivatives
        % d2udxdt = compactdiff(dudx,dt,3);
        % d2vdydt = compactdiff(dvdy,dt,3);
        % d_udwdz_dx = compactdiff(udwdz,dx,1);
        % d_vdwdz_dy = compactdiff(vdwdz,dy,2);
        % d2uvdydx = compactdiff(duvdy,dx,1);
        % d2uvdxdy = compactdiff(duvdx,dy,2);

        d_dwdz_dx = compactdiff(dwdz,dx,1);
        d_dwdz_dy = compactdiff(dwdz,dy,2);
    elseif strcmp(DIFFSCHEME,'soc')
        % %use new values for 1st derivatives to calculate higher time derivatives
        % d2udxdt    = socdiff(dudx,dt,3);
        % d2vdydt    = socdiff(dvdy,dt,3);
        % d_udwdz_dx = socdiff(udwdz,dx,1);
        % d_vdwdz_dy = socdiff(vdwdz,dy,2);
        % d2uvdydx   = socdiff(duvdy,dx,1);
        % d2uvdxdy   = socdiff(duvdx,dy,2);

        d_dwdz_dx = socdiff(dwdz,dx,1);
        d_dwdz_dy = socdiff(dwdz,dy,2);
    end
end
%% correct BC for higher order derivatives
%currently, none are used for the line integral formulations

%do higher order nested derivatives based on above BC derivatives
%calculate wall direction, normals, and derivatives at BC pts
%still going to have problems in 1-node acute points
%will crash if narrow channels close to domain wall face outwards (ie
%must be more than 2 points existing in direction BC points towards)
if strcmp(solver,'conservative2')
    for t=1:NT
    for n=1:length(bclist{t})
        flat = 1e-6;
        %disp(' xi  yi       s1       s2    nhat1    nhat2  N/S  E/W')

        %first point in list
        xi=bclist{t}{n}(1,1);
        yi=bclist{t}{n}(1,2);
        s = bclist{t}{n}(2,:) - bclist{t}{n}(end,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(1,:) = nhat;  %store unit normal for plotting
        %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)

        if nhat(2) < -flat      % north
            %use forward differences pointing downward
            %fprintf('     N')
            d2vdxdy   (xi, yi,:) = (3*  dvdx(xi, yi,:) -4*  dvdx(xi, yi-1,:) +   dvdx(xi, yi-2,:)) / (2*dy);
            %             d3vdx2dy  (xi, yi,:) = (3*d2vdx2(xi, yi,:) -4*d2vdx2(xi, yi-1,:) + d2vdx2(xi, yi-2,:)) / (2*dy);
            %             d3vdy2dy  (xi, yi,:) = (3*d2vdy2(xi, yi,:) -4*d2vdy2(xi, yi-1,:) + d2vdy2(xi, yi-2,:)) / (2*dy);
            d_dwdz_dy (xi, yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi, yi-1,:) +   dwdz(xi, yi-2,:)) / (2*dy);
            %             d_vdwdz_dy(xi, yi,:) = (3* vdwdz(xi, yi,:) -4* vdwdz(xi, yi-1,:) +  vdwdz(xi, yi-2,:)) / (2*dy);
            %             d2uvdxdy  (xi, yi,:) = (3* duvdx(xi, yi,:) -4* duvdx(xi, yi-1,:) +  duvdx(xi, yi-2,:)) / (2*dy);

        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
            %fprintf('     S')
            d2vdxdy   (xi, yi,:) = (-  dvdx(xi, yi+2,:) +4*  dvdx(xi, yi+1,:) -3*  dvdx(xi, yi,:)) / (2*dy);
            %             d3vdx2dy  (xi, yi,:) = (-d2vdx2(xi, yi+2,:) +4*d2vdx2(xi, yi+1,:) -3*d2vdx2(xi, yi,:)) / (2*dy);
            %             d3vdy2dy  (xi, yi,:) = (-d2vdy2(xi, yi+2,:) +4*d2vdy2(xi, yi+1,:) -3*d2vdy2(xi, yi,:)) / (2*dy);
            d_dwdz_dy (xi, yi,:) = (-  dwdz(xi, yi+2,:) +4*  dwdz(xi, yi+1,:) -3*  dwdz(xi, yi,:)) / (2*dy);
            %             d_vdwdz_dy(xi, yi,:) = (- vdwdz(xi, yi+2,:) +4* vdwdz(xi, yi+1,:) -3* vdwdz(xi, yi,:)) / (2*dy);
            %             d2uvdxdy  (xi, yi,:) = (- duvdx(xi, yi+2,:) +4* duvdx(xi, yi+1,:) -3* duvdx(xi, yi,:)) / (2*dy);

        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
            %fprintf('     0')
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
            %fprintf('     W\n')
            d2udydx   (xi ,yi,:) = (-  dudy(xi+2, yi,:) +4*  dudy(xi+1, yi,:) -3*  dudy(xi, yi,:)) / (2*dx);
            %             d3udx2dx  (xi ,yi,:) = (-d2udx2(xi+2, yi,:) +4*d2udx2(xi+1, yi,:) -3*d2udx2(xi, yi,:)) / (2*dx);
            %             d3udy2dx  (xi ,yi,:) = (-d2udy2(xi+2, yi,:) +4*d2udy2(xi+1, yi,:) -3*d2udy2(xi, yi,:)) / (2*dx);
            d_dwdz_dx (xi ,yi,:) = (-  dwdz(xi+2, yi,:) +4*  dwdz(xi+1, yi,:) -3*  dwdz(xi, yi,:)) / (2*dx);
            %             d_udwdz_dx(xi ,yi,:) = (- udwdz(xi+2, yi,:) +4* udwdz(xi+1, yi,:) -3* udwdz(xi, yi,:)) / (2*dx);
            %             d2uvdydx  (xi ,yi,:) = (- duvdy(xi+2, yi,:) +4* duvdy(xi+1, yi,:) -3* duvdy(xi, yi,:)) / (2*dx);

        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
            %fprintf('     E\n')
            d2udydx   (xi, yi,:) = (3*  dudy(xi, yi,:) -4*  dudy(xi-1, yi,:) +   dudy(xi-2, yi,:)) / (2*dx);
            %             d3udx2dx  (xi, yi,:) = (3*d2udx2(xi, yi,:) -4*d2udx2(xi-1, yi,:) + d2udx2(xi-2, yi,:)) / (2*dx);
            %             d3udy2dx  (xi, yi,:) = (3*d2udy2(xi, yi,:) -4*d2udy2(xi-1, yi,:) + d2udy2(xi-2, yi,:)) / (2*dx);
            d_dwdz_dx (xi ,yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi-1, yi,:) +   dwdz(xi-2, yi,:)) / (2*dx);
            %             d_udwdz_dx(xi ,yi,:) = (3* udwdz(xi, yi,:) -4* udwdz(xi-1, yi,:) +  udwdz(xi-2, yi,:)) / (2*dx);
            %             d2uvdydx  (xi ,yi,:) = (3* duvdy(xi, yi,:) -4* duvdy(xi-1, yi,:) +  duvdy(xi-2, yi,:)) / (2*dx);

        else
            %do nothing to x-derivatives
            %use central difference centered on this point
            %fprintf('     0\n')
        end
        %pause

        %middle pts
        for i=2:size(bclist{t}{n},1)-1
            xi=bclist{t}{n}(i,1);
            yi=bclist{t}{n}(i,2);
            s = bclist{t}{n}(i+1,:) - bclist{t}{n}(i-1,:);
            s = s/sqrt(s(1).^2+s(2).^2);
            nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
            nmlist{n}(i,:) = nhat;  %store unit normal for plotting
            %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)

            if nhat(2) < -flat      % north
                %use forward differences pointing downward
                %fprintf('     N')
                d2vdxdy   (xi, yi,:) = (3*  dvdx(xi, yi,:) -4*  dvdx(xi, yi-1,:) +   dvdx(xi, yi-2,:)) / (2*dy);
                %             d3vdx2dy  (xi, yi,:) = (3*d2vdx2(xi, yi,:) -4*d2vdx2(xi, yi-1,:) + d2vdx2(xi, yi-2,:)) / (2*dy);
                %             d3vdy2dy  (xi, yi,:) = (3*d2vdy2(xi, yi,:) -4*d2vdy2(xi, yi-1,:) + d2vdy2(xi, yi-2,:)) / (2*dy);
                d_dwdz_dy (xi, yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi, yi-1,:) +   dwdz(xi, yi-2,:)) / (2*dy);
                %             d_vdwdz_dy(xi, yi,:) = (3* vdwdz(xi, yi,:) -4* vdwdz(xi, yi-1,:) +  vdwdz(xi, yi-2,:)) / (2*dy);
                %             d2uvdxdy  (xi, yi,:) = (3* duvdx(xi, yi,:) -4* duvdx(xi, yi-1,:) +  duvdx(xi, yi-2,:)) / (2*dy);

            elseif nhat(2) > flat   % south
                %use forward differences pointing upward
                %fprintf('     S')
                d2vdxdy   (xi, yi,:) = (-  dvdx(xi, yi+2,:) +4*  dvdx(xi, yi+1,:) -3*  dvdx(xi, yi,:)) / (2*dy);
                %             d3vdx2dy  (xi, yi,:) = (-d2vdx2(xi, yi+2,:) +4*d2vdx2(xi, yi+1,:) -3*d2vdx2(xi, yi,:)) / (2*dy);
                %             d3vdy2dy  (xi, yi,:) = (-d2vdy2(xi, yi+2,:) +4*d2vdy2(xi, yi+1,:) -3*d2vdy2(xi, yi,:)) / (2*dy);
                d_dwdz_dy (xi, yi,:) = (-  dwdz(xi, yi+2,:) +4*  dwdz(xi, yi+1,:) -3*  dwdz(xi, yi,:)) / (2*dy);
                %             d_vdwdz_dy(xi, yi,:) = (- vdwdz(xi, yi+2,:) +4* vdwdz(xi, yi+1,:) -3* vdwdz(xi, yi,:)) / (2*dy);
                %             d2uvdxdy  (xi, yi,:) = (- duvdx(xi, yi+2,:) +4* duvdx(xi, yi+1,:) -3* duvdx(xi, yi,:)) / (2*dy);

            else                    % flat
                %do nothing to y-derivatives
                %use central difference centered on this point
                %fprintf('     0')
            end

            if nhat(1) > flat       %west
                %use forward differences pointing right
                %fprintf('     W\n')
                d2udydx   (xi ,yi,:) = (-  dudy(xi+2, yi,:) +4*  dudy(xi+1, yi,:) -3*  dudy(xi, yi,:)) / (2*dx);
                %             d3udx2dx  (xi ,yi,:) = (-d2udx2(xi+2, yi,:) +4*d2udx2(xi+1, yi,:) -3*d2udx2(xi, yi,:)) / (2*dx);
                %             d3udy2dx  (xi ,yi,:) = (-d2udy2(xi+2, yi,:) +4*d2udy2(xi+1, yi,:) -3*d2udy2(xi, yi,:)) / (2*dx);
                d_dwdz_dx (xi ,yi,:) = (-  dwdz(xi+2, yi,:) +4*  dwdz(xi+1, yi,:) -3*  dwdz(xi, yi,:)) / (2*dx);
                %             d_udwdz_dx(xi ,yi,:) = (- udwdz(xi+2, yi,:) +4* udwdz(xi+1, yi,:) -3* udwdz(xi, yi,:)) / (2*dx);
                %             d2uvdydx  (xi ,yi,:) = (- duvdy(xi+2, yi,:) +4* duvdy(xi+1, yi,:) -3* duvdy(xi, yi,:)) / (2*dx);

            elseif nhat(1) < -flat  %east
                %use forward differences pointing left
                %fprintf('     E\n')
                d2udydx   (xi, yi,:) = (3*  dudy(xi, yi,:) -4*  dudy(xi-1, yi,:) +   dudy(xi-2, yi,:)) / (2*dx);
                %             d3udx2dx  (xi, yi,:) = (3*d2udx2(xi, yi,:) -4*d2udx2(xi-1, yi,:) + d2udx2(xi-2, yi,:)) / (2*dx);
                %             d3udy2dx  (xi, yi,:) = (3*d2udy2(xi, yi,:) -4*d2udy2(xi-1, yi,:) + d2udy2(xi-2, yi,:)) / (2*dx);
                d_dwdz_dx (xi ,yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi-1, yi,:) +   dwdz(xi-2, yi,:)) / (2*dx);
                %             d_udwdz_dx(xi ,yi,:) = (3* udwdz(xi, yi,:) -4* udwdz(xi-1, yi,:) +  udwdz(xi-2, yi,:)) / (2*dx);
                %             d2uvdydx  (xi ,yi,:) = (3* duvdy(xi, yi,:) -4* duvdy(xi-1, yi,:) +  duvdy(xi-2, yi,:)) / (2*dx);

            else
                %do nothing to x-derivatives
                %use central difference centered on this point
                %fprintf('     0\n')
            end
            %pause
        end %middle points

        %last point in list
        xi=bclist{t}{n}(end,1);
        yi=bclist{t}{n}(end,2);
        s = bclist{t}{n}(1,:) - bclist{t}{n}(end-1,:);
        s = s/sqrt(s(1).^2+s(2).^2);
        nhat = [s(2),-s(1)];    %equivalent to dot(s,[0 0 1])
        nmlist{n}(end+1,:) = nhat;  %store unit normal for plotting
        %fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g',xi, yi, s, nhat)

        if nhat(2) < -flat      % north
            %use forward differences pointing downward
            %fprintf('     N')
            d2vdxdy   (xi, yi,:) = (3*  dvdx(xi, yi,:) -4*  dvdx(xi, yi-1,:) +   dvdx(xi, yi-2,:)) / (2*dy);
            %             d3vdx2dy  (xi, yi,:) = (3*d2vdx2(xi, yi,:) -4*d2vdx2(xi, yi-1,:) + d2vdx2(xi, yi-2,:)) / (2*dy);
            %             d3vdy2dy  (xi, yi,:) = (3*d2vdy2(xi, yi,:) -4*d2vdy2(xi, yi-1,:) + d2vdy2(xi, yi-2,:)) / (2*dy);
            d_dwdz_dy (xi, yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi, yi-1,:) +   dwdz(xi, yi-2,:)) / (2*dy);
            %             d_vdwdz_dy(xi, yi,:) = (3* vdwdz(xi, yi,:) -4* vdwdz(xi, yi-1,:) +  vdwdz(xi, yi-2,:)) / (2*dy);
            %             d2uvdxdy  (xi, yi,:) = (3* duvdx(xi, yi,:) -4* duvdx(xi, yi-1,:) +  duvdx(xi, yi-2,:)) / (2*dy);

        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
            %fprintf('     S')
            d2vdxdy   (xi, yi,:) = (-  dvdx(xi, yi+2,:) +4*  dvdx(xi, yi+1,:) -3*  dvdx(xi, yi,:)) / (2*dy);
            %             d3vdx2dy  (xi, yi,:) = (-d2vdx2(xi, yi+2,:) +4*d2vdx2(xi, yi+1,:) -3*d2vdx2(xi, yi,:)) / (2*dy);
            %             d3vdy2dy  (xi, yi,:) = (-d2vdy2(xi, yi+2,:) +4*d2vdy2(xi, yi+1,:) -3*d2vdy2(xi, yi,:)) / (2*dy);
            d_dwdz_dy (xi, yi,:) = (-  dwdz(xi, yi+2,:) +4*  dwdz(xi, yi+1,:) -3*  dwdz(xi, yi,:)) / (2*dy);
            %             d_vdwdz_dy(xi, yi,:) = (- vdwdz(xi, yi+2,:) +4* vdwdz(xi, yi+1,:) -3* vdwdz(xi, yi,:)) / (2*dy);
            %             d2uvdxdy  (xi, yi,:) = (- duvdx(xi, yi+2,:) +4* duvdx(xi, yi+1,:) -3* duvdx(xi, yi,:)) / (2*dy);

        else                    % flat
            %do nothing to y-derivatives
            %use central difference centered on this point
            %fprintf('     0')
        end

        if nhat(1) > flat       %west
            %use forward differences pointing right
            %fprintf('     W\n')
            d2udydx   (xi ,yi,:) = (-  dudy(xi+2, yi,:) +4*  dudy(xi+1, yi,:) -3*  dudy(xi, yi,:)) / (2*dx);
            %             d3udx2dx  (xi ,yi,:) = (-d2udx2(xi+2, yi,:) +4*d2udx2(xi+1, yi,:) -3*d2udx2(xi, yi,:)) / (2*dx);
            %             d3udy2dx  (xi ,yi,:) = (-d2udy2(xi+2, yi,:) +4*d2udy2(xi+1, yi,:) -3*d2udy2(xi, yi,:)) / (2*dx);
            d_dwdz_dx (xi ,yi,:) = (-  dwdz(xi+2, yi,:) +4*  dwdz(xi+1, yi,:) -3*  dwdz(xi, yi,:)) / (2*dx);
            %             d_udwdz_dx(xi ,yi,:) = (- udwdz(xi+2, yi,:) +4* udwdz(xi+1, yi,:) -3* udwdz(xi, yi,:)) / (2*dx);
            %             d2uvdydx  (xi ,yi,:) = (- duvdy(xi+2, yi,:) +4* duvdy(xi+1, yi,:) -3* duvdy(xi, yi,:)) / (2*dx);

        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
            %fprintf('     E\n')
            d2udydx   (xi, yi,:) = (3*  dudy(xi, yi,:) -4*  dudy(xi-1, yi,:) +   dudy(xi-2, yi,:)) / (2*dx);
            %             d3udx2dx  (xi, yi,:) = (3*d2udx2(xi, yi,:) -4*d2udx2(xi-1, yi,:) + d2udx2(xi-2, yi,:)) / (2*dx);
            %             d3udy2dx  (xi, yi,:) = (3*d2udy2(xi, yi,:) -4*d2udy2(xi-1, yi,:) + d2udy2(xi-2, yi,:)) / (2*dx);
            d_dwdz_dx (xi ,yi,:) = (3*  dwdz(xi, yi,:) -4*  dwdz(xi-1, yi,:) +   dwdz(xi-2, yi,:)) / (2*dx);
            %             d_udwdz_dx(xi ,yi,:) = (3* udwdz(xi, yi,:) -4* udwdz(xi-1, yi,:) +  udwdz(xi-2, yi,:)) / (2*dx);
            %             d2uvdydx  (xi ,yi,:) = (3* duvdy(xi, yi,:) -4* duvdy(xi-1, yi,:) +  duvdy(xi-2, yi,:)) / (2*dx);

        else
            %do nothing to x-derivatives
            %use central difference centered on this point
            %fprintf('     0\n')
        end
        %pause
    end %end recalculation of higher derivatives on BC
    end %for k=1:NT
end %if conservative2
%keyboard

%% Calculate N-S field function        
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
    error('unknown solver');
end

%keyboard

%can probably prune the list of derivatives for line integrals - many terms
%are for the source term only
clear dudt dvdt dudx dvdx dudy dvdy d2udx2 d2udy2 d2vdx2 d2vdy2
clear dfdx dgdy uu uv vv udwdz vdwdz
clear d2udydx d3udx2dx d3udy2dx d2vdxdy d3vdx2dy d3vdy2dy d2udxdt d2vdydt
clear d2uudx2 d2uvdydx d_udwdz_dx d2uvdxdy d2vvdy2 d_vdwdz_dy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     disp(' xi  yi       s1       s2    nhat1    nhat2       Ap       Ae       Aw       An       As        f        g        S')


ktot = 0;
relres(ktot+1) = tol;   %change in pressure frame to frame
relcon(ktot+1) = 1;     %change in pressure residual frame to frame
relresf(ktot+1) = 1;    %error between dP/dx and F
relresg(ktot+1) = 1;    %error between dP/dy and G

%figure(942),hold off,semilogy(0,0,'w')

pause(0.001)

%solve for pressures
fprintf('     t       k        resid  convergence  no-visit  avg-visit        time\n')
fprintf('  ----  ------  -----------  -----------  --------  ---------  ----------\n')
% fprintf('  %4i  %6i  %11.4g  %4i  %10.3f\n',0,(iter(1)-1)*restarts+iter(2),relres,flag,toc)

%keyboard

%% setup paths
    %diagonal of domain
    DD = sqrt( (NI-1)^2 + (NJ-1)^2 );

    %center of domain
    cX = (1+NI)/2;
    cY = (1+NJ)/2;
    
    RR = DD/2 * RSIZE; %RSIZE ~1.5, >1 
  
    %for evenly spaced nodes and paths:
    theta = linspace(0,2*pi,NR+1);          %don't need repeated endpoint at 2pi
    boundX = round(RR*cos(theta(1:end-1))+cX);  %round to integer grid locations
    boundY = round(RR*sin(theta(1:end-1))+cY);

%     %for random paths:
%     theta = 2*pi*rand(NR,1);
%     boundX = round(RR*cos(theta(1:end))+cX);  %round to integer grid locations
%     boundY = round(RR*sin(theta(1:end))+cY);

%% integrations
for t=1:NT


    
    pold   = zeros(NI,NJ);  %at first iteration will only contain bc values
    
    %integrate along BC's - begin each BC with P=0
    for n=1:length(bclist{t})
        dxlist = dx*[ 0 ; bclist{t}{n}(2:end,1) - bclist{t}{n}(1:end-1,1) ];
        dylist = dy*[ 0 ; bclist{t}{n}(2:end,2) - bclist{t}{n}(1:end-1,2) ];
        dnlist = (dxlist.^2 + dylist.^2).^2;
        thlist = atan2( dylist, dxlist );
        
        
        
        %k=1
%             kx = bclist{t}{n}(1,1);
%             ky = bclist{t}{n}(1,2);
            ptotal = pold( bclist{t}{n}(1,1), bclist{t}{n}(1,2) );
        for k=1:length(bclist{t}{n})
%             kx = bclist{t}{n}(k,1);
%             ky = bclist{t}{n}(k,2);
            %calculate new pressure at k
%             ptotal = ptotal + ( cos(thlist(k))*( f( bclist{t}{n}(k,1), bclist{t}{n}(k,2), t ) + f( bclist{t}{n}(k-1,1), bclist{t}{n}(k-1,2) ,t ) )/2 ...
%                                +sin(thlist(k))*( g( bclist{t}{n}(k,1), bclist{t}{n}(k,2), t ) + g( bclist{t}{n}(k-1,1), bclist{t}{n}(k-1,2) ,t ) )/2 )*dnlist(k);
            %rewrite duplicate values
            pold(bclist{t}{n}(1,k), bclist{t}{n}(2,k))   = ptotal;
            %pcount(bclist{t}{n}(k,1), bclist{t}{n}(k,2)) = 1;
        end
    end
    
    %iterate this section till convergence
    
    for q=1:QMAX
        ktot = ktot +1;
        pnew   = zeros(NI,NJ);
        pcount = zeros(NI,NJ);
        

        %for m=1:NR-1
        for m=1:NR
            for n=[1:m-1,m+1:NR]
                %for n=m+1:NR
                %setup path
                pathDX = abs(boundX(n)-boundX(m))+1;   %number of grid pts between start and end
                pathDY = abs(boundY(n)-boundY(m))+1;
                pathDL = max(pathDX,pathDY);
                %find grid points along path
                % pathX = round(linspace(boundX(m),boundX(n),pathDL));
                % pathY = round(linspace(boundY(m),boundY(n),pathDL));
                % pathX = round(lspace(boundX(m),boundX(n),pathDL));
                % pathY = round(lspace(boundY(m),boundY(n),pathDL));
                pathX = round( [boundX(m)+(0:pathDL-2)*(boundX(n)-boundX(m))/(pathDL-1) boundX(n)] );
                pathY = round( [boundY(m)+(0:pathDL-2)*(boundY(n)-boundY(m))/(pathDL-1) boundY(n)] );
                                
                %constrain pts to within borders of domain
                %             temp = find(pathX >= 1 & pathX <= NI);
                %             pathX = pathX(temp);
                %             pathY = pathY(temp);
                %             temp = find(pathY >= 1 & pathY <= NJ);
                %             pathX = pathX(temp);
                %             pathY = pathY(temp);
                temp = find(pathX >= 1 & pathX <= NI & pathY >= 1 & pathY <= NJ);
                pathX = pathX(temp);
                pathY = pathY(temp);
                %determine which pts are within BC's
                %maskPts = bcmask(sub2ind(size(bcmask),pathX,pathY,t));	%extract bc info from pts along path
                %maskPts = bcmask(s2i3([NI,NJ,NT],pathX,pathY,t));	%extract bc info from pts along path
                maskPts = bcmask(NI*NJ*(t-1)+NI*(pathY-1)+pathX);	%extract bc info from pts along path
                %             intPts  = [0, maskPts(2:end) & maskPts(1:end-1)];       %determine which have previous valid pts
                intPts  = ( maskPts(2:end) & maskPts(1:end-1) );       %determine which have previous valid pts
                %calculate distance between pts, and angle of step
                %             dxlist = dx*[ 0 , pathX(2:end) - pathX(1:end-1) ];
                %             dylist = dy*[ 0 , pathY(2:end) - pathY(1:end-1) ];
                dxlist = dx*( pathX(2:end) - pathX(1:end-1) );
                dylist = dy*( pathY(2:end) - pathY(1:end-1) );
                dnlist = sqrt(dxlist.^2 + dylist.^2);
                thlist = atan2( dylist, dxlist );

                %ptotal = 0;
                if numel(pathX)~=0
                    k=1;
                    ptotal = pold( pathX(k), pathY(k));                         %reset line total to previous value (probably 0)
                    %pnew(pathX(k), pathY(k)) = pnew(pathX(k), pathY(k));        %initialize ptotal
                    pnew(pathX(k), pathY(k)) = pnew(pathX(k), pathY(k)) + ptotal;        %initialize ptotal
                    pcount(pathX(k), pathY(k)) = pcount(pathX(k), pathY(k)) + 1;    %increment visit counter
                end
                for k=2:length(pathX)   %integrate other positions
                    if intPts(k-1) %need to integrate to next point
                        ptotal = ptotal + ( cos(thlist(k-1))*( f( pathX(k), pathY(k), t ) + f( pathX(k-1), pathY(k-1) ,t ) )/2 ...
                                           +sin(thlist(k-1))*( g( pathX(k), pathY(k), t ) + g( pathX(k-1), pathY(k-1) ,t ) )/2 )*dnlist(k-1);
                        pnew(pathX(k),pathY(k)) = pnew(pathX(k),pathY(k)) + ptotal;
                        pcount(pathX(k),pathY(k)) = pcount(pathX(k),pathY(k)) + 1;
                    else %either bc or outside valid region
                        ptotal = pold( pathX(k), pathY(k));                         %reset line total to previous value (probably 0)
                        pnew(pathX(k), pathY(k)) = pnew(pathX(k), pathY(k)) + ptotal;   %add ptotal to field
                        pcount(pathX(k), pathY(k)) = pcount(pathX(k), pathY(k)) + 1;    %increment visit counter
                    end

                end %path

                %             ptotal = 0;
                %             for k=1:length(pathX)   %integrate other positions
                %                 if intPts(k) %need to integrate to next point
                %                 	ptotal = ptotal + ( cos(thlist(k))*( f( pathX(k), pathY(k), t ) + f( pathX(k-1), pathY(k-1) ,t ) )/2 ...
                %                                        +sin(thlist(k))*( g( pathX(k), pathY(k), t ) + g( pathX(k-1), pathY(k-1) ,t ) )/2 )*dnlist(k);
                %                     pnew(pathX(k),pathY(k)) = pnew(pathX(k),pathY(k)) + ptotal;
                %                     pcount(pathX(k),pathY(k)) = pcount(pathX(k),pathY(k)) + 1;
                %                 else %either bc or outside valid region
                %                 	ptotal = pold( pathX(k), pathY(k));                         %reset line total to previous value (probably 0)
                %                     pnew(pathX(k), pathY(k)) = pnew(pathX(k), pathY(k)) + ptotal;   %add ptotal to field
                %                     pcount(pathX(k), pathY(k)) = pcount(pathX(k), pathY(k)) + 1;    %increment visit counter
                %                 end
                %
                %             end %path

            end %all possible end pts from starting point
        end %all possible starting points

        %     figure(801),imagesc(pcount.'),axis image xy,title('visits'),colorbar

        %make sure every cell is visited at least once
        temp = find(pcount==0);
        no_visits  = length(temp);
        avg_visits = mean(mean(pcount));
%         disp(['The number of unvisited cells is ' num2str(no_visits) ]);
%         disp(['The average visit count is ' num2str(avg_visits) ])
        pcount(temp) = 1;

         %pcount( find(pcount==0) ) = 1;
        pnew = pnew./pcount;    %take average of all visits
        
        relres(ktot+1) = sqrt(sum( (pnew(:)-pold(:)).^2 )) / sqrt(sum( pnew(:).^2 ));
        relcon(ktot+1) = relres(ktot+1)/relres(ktot);
        
        dpdx = socdiff_bc(pnew,dx,1,bcmask(:,:,t));
        dpdy = socdiff_bc(pnew,dy,2,bcmask(:,:,t));
        
        resf = dpdx-f(:,:,t);
        resg = dpdy-g(:,:,t);
        
        relresf(ktot+1) = sqrt(sum( resf(:).^2 ) ) ./ sqrt(sum( dpdx(:).^2));
        relresg(ktot+1) = sqrt(sum( resg(:).^2 ) ) ./ sqrt(sum( dpdy(:).^2));

        figure(801),semilogy(1:ktot,relres(2:ktot+1),'o-',...
                             1:ktot,1-relcon(2:ktot+1),'o-',...
                             1:ktot,relresf(2:ktot+1),'o-',...
                             1:ktot,relresg(2:ktot+1),'o-'),legend('residual','1-convergence')

        flag=0;iter=[1,1];resvec=0;%relres=0;
        %convergence information
        fprintf('  %4i  %6i  %11.4g  %11.4g  %8g  %9g  %10.3f\n',t,q,relres(ktot+1),relcon(ktot+1),no_visits,avg_visits,toc)
        %keyboard
        
        pmin = min(min(pnew));
        pmax = max(max(pnew));

            figure(802),imagesc(1100/133.3*pold.',[1100/133.3*pmin 1100/133.3*pmax]),axis image xy,title('old pressure field'),colorbar;
            figure(803),imagesc(1100/133.3*pnew.',[1100/133.3*pmin 1100/133.3*pmax]),axis image xy,title('new pressure field'),colorbar;
            figure(804),imagesc(((pnew-pold)./(pmax-pmin)).'),axis image xy,title('relative error'),colorbar;
%             figure(805),imagesc(resf.'),axis image xy,title('residual on dP/dx'),colorbar;
%             figure(806),imagesc(resg.'),axis image xy,title('residual on dP/dy'),colorbar;
%         pause(0.01)
       

        %pold = pnew;    %copy pnew over pold in prep for next iteration
        AA = 1.5;
        pold = pold + AA*(pnew-pold);    %copy pnew over pold in prep for next iteration
        
        %finish iteration section
        if relres(ktot+1)<tol && relcon(ktot+1)>tolcon && relcon(ktot+1)<1
            break
        end

    end %iterations around bc
    
    %store final solution into pressure table
    p(:,:,t) = pold;    %pold should hold final iteration
    
    %use an offset to peg the reference pressure over time... 
    %is this okay to do?
    if nargin>=10
        %pause
        deltaP = Pref(t) - p(Px,Py,t);
        p(:,:,t) = p(:,:,t) + deltaP * ones(NI,NJ);
    end

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
       
     %ktot = ktot + length(resvec)-1;
    

end %proceed to next time step
% catch
%     disp(lasterr)
%     keyboard
% end

warning('on','MATLAB:Axes:NegativeDataInLogAxis');


function N=s2i(siz,I,J)
%faster version of sub2ind
%converts I,J subscript into linear index into array
N=siz(1)*(J-1)+I;

function N=s2i3(siz,I,J,K)
%faster version of sub2ind
%converts I,J,K subscript into linear index into array
N=siz(2)*siz(1)*(K-1)+siz(1)*(J-1)+I;

function y = lspace(d1, d2, n)
%LINSPACE Linearly spaced vector.
%   LINSPACE(X1, X2) generates a row vector of 100 linearly
%   equally spaced points between X1 and X2.
%
%   LINSPACE(X1, X2, N) generates N points between X1 and X2.
%   For N < 2, LINSPACE returns X2.
%
%   Class support for inputs X1,X2:
%      float: double, single
%
%   See also LOGSPACE, :.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.12.4.1 $  $Date: 2004/07/05 17:01:20 $
%  slightly adapted for speed by John Charonko 2007/11/29

y = [d1+(0:n-2)*(d2-d1)/(n-1) d2];


