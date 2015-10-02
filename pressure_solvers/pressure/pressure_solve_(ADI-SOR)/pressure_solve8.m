function [p]=pressure_solve(x,y,u,v,Re,dt,flowtype,bclist,Pref,Px,Py)
%function [S]=pressure_solve(x,y,u,v,Re,dt,[flowtype],[bclist],[Pref, Px, Py])
%
% Returns p(x,y,t) from u(x,y,t) and v(x,y,t) using the ADI iterative 
% method with SOR to solve the pressure Poisson equation:
%   del^2(P) = S(x,u)
%     S(x,y) = df/dx + dg/dy
%   From the N-S equations:
%     f = 1/Re * (d2u/dx2+d2u/dy2) - ( du/dt + u*du/dx + v*du/dy )
%     g = 1/Re * (d2v/dx2+d2v/dy2) - ( dv/dt + u*dv/dx + v*dv/dy )
%   but for incompressible flow, the source term simplifies to:
%      S(x,y) = du_i/dx_j*du_j/dx_i
%
% FLOWTYPE must be either "viscous" or "ideal"
% The solver leaves off the viscous terms for "ideal" flow solutions.
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

% updated derivative calc on boundaries of domains...
% updated derivative calc on given boundary conditions
% added pressure Dirichlet BC Pref in interior at Px and Py

% based on tafti-hw9_1.m
% solve the Poisson equation del^2(phi)=S(x,y)
% S(x,y) = -4*a*b(1-b*(x^2+y^2))*exp(-b*(x^2+y^2)), with a=b=1;
% over 0<=x<=1, 0<=y<=1
% using the ADI iterative method with SOR

% try

figure(1234),quiver(x,y,u(:,:,1),v(:,:,1)),pause

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
            disp('PREF, PX, and PY are defined')
            Pref = Pref(:); %trick to make Pref vertical
            if (length(Pref)==1)
                Pref = Pref*ones(length(Pref),1);
            else
                if length(Pref)~=size(u,3)
                    error('PREF must have length 1 or equal to number of flow fields in U and V');
                end
            end
        end
        disp('BCLIST is defined')
        if iscell(bclist)
            disp('BCLIST is a cell')
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
    flowtype = 'viscous'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
tic
max_iterations  = 100000;
loop_iterations = 100;
wglobal = 0.05;

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
An = zeros(NI,NJ);  
As = zeros(NI,NJ);
Aw = zeros(NI,NJ);
Ae = zeros(NI,NJ);
Ap = zeros(NI,NJ);
%store pressure in phi during solution
phi = zeros(NI,NJ,3);
%source term plus LHS terms not used in tridiagonal solution
RHSi = zeros(NI,1);
RHSj = zeros(NJ,1);
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
dudt(:,:,1    ) = (u(:,:,2  ) - u(:,:,1    )) /    dt ;
dudt(:,:,2:T-1) = (u(:,:,3:T) - u(:,:,1:T-2)) / (2*dt);
dudt(:,:,    T) = (u(:,:,  T) - u(:,:,  T-1)) /    dt ;

dvdt(:,:,1    ) = (v(:,:,2  ) - v(:,:,1    )) /    dt ;
dvdt(:,:,2:T-1) = (v(:,:,3:T) - v(:,:,1:T-2)) / (2*dt);
dvdt(:,:,    T) = (v(:,:,  T) - v(:,:,  T-1)) /    dt ;

dudx(1    ,:,:) = ( -u(3  ,:,:) +4*u( 2 ,:,:) -3*u(1    ,:,:)) / (2*dx);
dudx(2:M-1,:,:) = (  u(3:M,:,:)               -  u(1:M-2,:,:)) / (2*dx);
dudx(M    ,:,:) = (3*u(M  ,:,:) -4*u(M-1,:,:) +  u(  M-2,:,:)) / (2*dx);

dvdx(1    ,:,:) = ( -v(3  ,:,:) +4*v( 2 ,:,:) -3*v(1    ,:,:)) / (2*dx);
dvdx(2:M-1,:,:) = (  v(3:M,:,:)               -  v(1:M-2,:,:)) / (2*dx);
dvdx(M    ,:,:) = (3*v(M  ,:,:) -4*v(M-1,:,:) +  v(  M-2,:,:)) / (2*dx);

dudy(:,1    ,:) = ( -u(:,3  ,:) +4*u(:, 2 ,:) -3*u(:,1    ,:)) / (2*dx);
dudy(:,2:N-1,:) = (  u(:,3:N,:)               -  u(:,1:N-2,:)) / (2*dy);
dudy(:,    N,:) = (3*u(:,  N,:) -4*u(:,N-1,:) +  u(:,  N-2,:)) / (2*dy);

dvdy(:,1    ,:) = ( -v(:,3  ,:) +4*v(:, 2 ,:) -3*v(:,1    ,:)) / (2*dx);
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

%calculate wall direction, normals, and derivatives at BC pts
%still going to have problems in 1-node acute points
%will crash if narrow channels close to domain wall face outwards (ie
%must be more than 2 points existing in direction BC points towards)
for n=1:length(bclist)
    flat = 1e-6;
    %disp(' xi  yi       s1       s2    nhat1    nhat2  N/S  E/W')
    
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
            %fprintf('     N')
            dudy  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi, yi-1,:) + u(xi, yi-2,:)) / (2*dy);
            dvdy  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi, yi-1,:) + v(xi, yi-2,:)) / (2*dy);
            d2udy2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi, yi-1,:) + u(xi, yi-2,:)) / (dy^2);
            d2vdy2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi, yi-1,:) + v(xi, yi-2,:)) / (dy^2);
            
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
            %fprintf('     S')
            dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dx);
            dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dx);
            d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
            d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
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

        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
            %fprintf('     E\n')
            dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
            dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
            d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
            d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
            %fprintf('     0\n')
        end
        %pause
        
    %middle pts
    for i=2:size(bclist{n},1)-1
        xi=bclist{n}(i,1);
        yi=bclist{n}(i,2);
        s = bclist{n}(i+1,:) - bclist{n}(i-1,:);
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
            
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
            %fprintf('     S')
            dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dx);
            dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dx);
            d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
            d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
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

        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
            %fprintf('     E\n')
            dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
            dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
            d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
            d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
            %fprintf('     0\n')
        end
        %pause
    end %middle points
        
    %last point in list
        xi=bclist{n}(end,1);
        yi=bclist{n}(end,2);
        s = bclist{n}(1,:) - bclist{n}(end-1,:);
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
            
        elseif nhat(2) > flat   % south
            %use forward differences pointing upward
            %fprintf('     S')
            dudy  (xi, yi,:) = (-u(xi, yi+2,:) +4*u(xi, yi+1,:) -3*u(xi, yi,:)) / (2*dx);
            dvdy  (xi, yi,:) = (-v(xi, yi+2,:) +4*v(xi, yi+1,:) -3*v(xi, yi,:)) / (2*dx);
            d2udy2(xi, yi,:) = ( u(xi, yi+2,:) -2*u(xi, yi+1,:) +  u(xi, yi,:)) / (dy^2);
            d2vdy2(xi, yi,:) = ( v(xi, yi+2,:) -2*v(xi, yi+1,:) +  v(xi, yi,:)) / (dy^2);
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

        elseif nhat(1) < -flat  %east
            %use forward differences pointing left
            %fprintf('     E\n')
            dudx  (xi, yi,:) = (3*u(xi, yi,:) -4*u(xi-1, yi,:) + u(xi-2, yi,:)) / (2*dx);
            dvdx  (xi, yi,:) = (3*v(xi, yi,:) -4*v(xi-1, yi,:) + v(xi-2, yi,:)) / (2*dx);
            d2udx2(xi, yi,:) = (  u(xi, yi,:) -2*u(xi-1, yi,:) + u(xi-2, yi,:)) / (dx^2);
            d2vdx2(xi, yi,:) = (  v(xi, yi,:) -2*v(xi-1, yi,:) + v(xi-2, yi,:)) / (dx^2);
        else            
            %do nothing to x-derivatives
            %use central difference centered on this point
            %fprintf('     0\n')
        end
        %pause
    figure(2345)
    plot(bclist{n}(:,1),bclist{n}(:,2))
    hold on
    quiver(bclist{n}(:,1),bclist{n}(:,2),nmlist{n}(:,1),nmlist{n}(:,2))
    hold off
    %keyboard
end %end recalculation of derivatives on BC



%Calculate pressure gradients        
if strcmp(flowtype,'ideal')
    f = - ( dudt + u.*dudx + v.*dudy );
    g = - ( dvdt + u.*dvdx + v.*dvdy );
elseif strcmp(flowtype,'viscous')
    f = 1/Re * (d2udx2+d2udy2) - ( dudt + u.*dudx + v.*dudy );
    g = 1/Re * (d2vdx2+d2vdy2) - ( dvdt + u.*dvdx + v.*dvdy );
else
    error('unknown flowtype');
end

%%%%%%%%%%
%Calculate source terms
%%%%%%%%%%
% dfdx(1    ,:,:) = (f(2  ,:,:) - f(1    ,:,:)) /    dx ;
% dfdx(2:M-1,:,:) = (f(3:M,:,:) - f(1:M-2,:,:)) / (2*dx);
% dfdx(M    ,:,:) = (f(M  ,:,:) - f(  M-1,:,:)) /    dx ;
% 
% dgdy(:,1    ,:) = (g(:,2  ,:) - g(:,1    ,:)) /    dy ;
% dgdy(:,2:N-1,:) = (g(:,3:N,:) - g(:,1:N-2,:)) / (2*dy);
% dgdy(:,    N,:) = (g(:,  N,:) - g(:,  N-1,:)) /    dy ;
% 
% Source = dfdx + dgdy;

%dwdz correction based on continuity
dwdz = - (dudx + dvdy);

% Source = -( dudx.*dudx + dudy.*dvdx ...
%            +dvdx.*dudy + dvdy.*dvdy ); 

Source = -( dudx.*dudx + dudy.*dvdx ...
           +dvdx.*dudy + dvdy.*dvdy ...
           +dwdz.*dwdz );        
       
%keyboard;
       
clear dudt dvdt dudx dvdx dudy dvdy d2udx2 d2udy2 d2vdx2 d2vdy2
clear dfdx dgdy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use SOC for middle of region
% the problem has been discritized so that
% Ap*phi_i,j    + 
%  Ae*phi_i+1,j + Aw*phi_i-1,j  + 
%  An*phi_i,j   + As*phi_i,j-1   =  S_i,j
% 

%default discretization
Ap(1:NI,1:NJ)   = -2*(1/dx^2+1/dy^2);
Ae(1:NI,1:NJ)   = 1/dx^2;
Aw(1:NI,1:NJ)   = 1/dx^2;
An(1:NI,1:NJ)   = 1/dy^2;
As(1:NI,1:NJ)   = 1/dy^2;
 S(1:NI,1:NJ,:) = Source(1:NI,1:NJ,:);


%calculate wall direction, normals, and BC
for n=1:length(bclist)
%     disp(' xi  yi       s1       s2    nhat1    nhat2       Ap       Ae       Aw       An       As        f        g        S')
    
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

%         fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',xi, yi, s, nhat, Ap(xi,yi), Ae(xi,yi) ,Aw(xi,yi) ,An(xi,yi) ,As(xi,yi), f(xi,yi,1),  g(xi,yi,1), S(xi,yi,1))
%         pause
        
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
%         fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',xi, yi, s, nhat, Ap(xi,yi), Ae(xi,yi) ,Aw(xi,yi) ,An(xi,yi) ,As(xi,yi), f(xi,yi,1),  g(xi,yi,1), S(xi,yi,1))
%         pause        
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
%         fprintf('%3i %3i %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',xi, yi, s, nhat, Ap(xi,yi), Ae(xi,yi) ,Aw(xi,yi) ,An(xi,yi) ,As(xi,yi), f(xi,yi,1),  g(xi,yi,1), S(xi,yi,1))
%         pause
end %for n=1:length(bclist) 

%Bad idea - convergence is terrible - just cheat and shift the 
if nargin>=9
    %Fixed pressure measurement at a point
    Ap(Px,Py)   = 1;
    Ae(Px,Py)   = 0;
    Aw(Px,Py)   = 0;
    An(Px,Py)   = 0;
    As(Px,Py)   = 0;
     S(Px,Py,:) = Pref;
end


%need this?
%solution copies previous solution to initial guess on each iteration, so
%we need to make sure there is something in 2 and 3
phi(:,:,2)=phi(:,:,3);  %fill in initial values for 

ktot = 0;
figure(942),hold off,semilogy(0,0,'w')
%solve for pressures

fprintf('     t       k        resid       lambda      w        time\n')
fprintf('  ----  ------  -----------  -----------  -----  ----------\n')

for t=1:NT

    %w=1.0;
    w=wglobal;  %defined at top
    k=1;
    resid=1;   %dummy value to start loop
    lambda=0;   %dummy value to start loop
    
    %find the rms residual and errors
    rho_tot = 0;
    e_tot = 0;
    rho_rms(k) = 1;
    e_rms(k)   = 1;

    resid  = rho_rms(k);

    while (k < max_iterations) & (resid > 1e-5) & (resid < 1e5) & (lambda < (1-1e-5))
        k=k+1;
        %advance the solution (erases previous timestep)
        phi(:,:,1) = phi(:,:,3);    %should preserve BC
        
        
        %iterate up j to calculate phi(i,j,k+1/2), take lines in i direction
          
        j=1;    %south
%        RHSi(2:NI-1) = (1-w)*Ap(2:NI-1,j).*phi(2:NI-1,j,1)-w*An(2:NI-1,j).*phi(2:NI-1,j+1,1)+w*g(2:NI-1,j,t);
%        RHSi(1 )     = (1-w)*Ap(1     ,j).*phi(1     ,j,1)-w*An(1     ,j).*phi(1     ,j+1,1)+w*f(1     ,j,t);   %west BC
%        RHSi(NI)     = (1-w)*Ap(NI    ,j).*phi(NI    ,j,1)-w*An(NI    ,j).*phi(NI    ,j+1,1)+w*g(NI    ,j,t);   %east BC
         RHSi(1:NI)   = (1-w)*Ap(1:NI  ,j).*phi(1:NI  ,j,1)-w*An(1:NI  ,j).*phi(1:NI  ,j+1,1)+w*S(1:NI  ,j,t);
        %solve for k+1/2
         phi(1:NI,j,2) = TDS1(1,NI,w*Aw(1:NI,j),Ap(1:NI,j),w*Ae(1:NI,j),RHSi(1:NI));        
        for j=2:NJ-1
%             RHSi(2:NI-1) = (1-w)*Ap(2:NI-1,j).*phi(2:NI-1,j,1)-w*An(2:NI-1,j).*phi(2:NI-1,j+1,1)-w*As(2:NI-1,j).*phi(2:NI-1,j-1,2)+w*S(2:NI-1,j,t);
%             RHSi(1 )     = (1-w)*Ap(1     ,j).*phi(1     ,j,1)-w*An(1     ,j).*phi(1     ,j+1,1)-w*As(1     ,j).*phi(1     ,j-1,2)+w*f(1     ,j,t);   %west BC
%             RHSi(NI)     = (1-w)*Ap(NI    ,j).*phi(NI    ,j,1)-w*An(NI    ,j).*phi(NI    ,j+1,1)-w*As(NI    ,j).*phi(NI    ,j-1,2)+w*f(NI    ,j,t);   %east BC
            RHSi(1:NI)   = (1-w)*Ap(1:NI  ,j).*phi(1:NI  ,j,1)-w*An(1:NI  ,j).*phi(1:NI  ,j+1,1)-w*As(1:NI  ,j).*phi(1:NI  ,j-1,2)+w*S(1:NI  ,j,t);
            %solve for k+1/2
            phi(1:NI,j,2) = TDS1(1,NI,w*Aw(1:NI,j),Ap(1:NI,j),w*Ae(1:NI,j),RHSi(1:NI));
        end
        j=NJ;   %north
%          RHSi(2:NI-1) = (1-w)*Ap(2:NI-1,j).*phi(2:NI-1,j,1)-w*As(2:NI-1,j).*phi(2:NI-1,j-1,2)+w*g(2:NI-1,j,t);
%          RHSi(1 )     = (1-w)*Ap(1     ,j).*phi(1     ,j,1)-w*As(1     ,j).*phi(1     ,j-1,2)+w*g(1     ,j,t);   %west BC
%          RHSi(NI)     = (1-w)*Ap(NI    ,j).*phi(NI    ,j,1)-w*As(NI    ,j).*phi(NI    ,j-1,2)+w*f(NI    ,j,t);   %east BC
         RHSi(1:NI)   = (1-w)*Ap(1:NI  ,j).*phi(1:NI  ,j,1)-w*As(1:NI  ,j).*phi(1:NI  ,j-1,1)+w*S(1:NI  ,j,t);
         %solve for k+1/2
         phi(1:NI,j,2) = TDS1(1,NI,w*Aw(1:NI,j),Ap(1:NI,j),w*Ae(1:NI,j),RHSi(1:NI));        
        
        %iterate across i to calculate phi(i,j,k+1/2), take lines in j direction
        i=1;    %west
%          RHSj(2:NJ-1) = (1-w)*Ap(i,2:NJ-1).*phi(i,2:NJ-1,2)-w*Ae(i,2:NJ-1).*phi(i+1,2:NJ-1,2)+w*f(i,2:NJ-1,t);
%          RHSj(1 )     = (1-w)*Ap(i,     1).*phi(i,     1,2)-w*Ae(i,     1).*phi(i+1,     1,2)+w*f(i,     1,t);   %south BC
%          RHSj(NJ)     = (1-w)*Ap(i,    NJ).*phi(i,    NJ,2)-w*Ae(i,    NJ).*phi(i+1,    NJ,2)+w*g(i,    NJ,t);   %north BC
         RHSj(1:NJ)   = (1-w)*Ap(i,  1:NJ).*phi(i,  1:NJ,2)-w*Ae(i,  1:NJ).*phi(i+1,  1:NJ,2)+w*S(i,  1:NJ,t);
         %solve for k+1
         phi(i,1:NJ,3) = TDS1(1,NJ,w*As(i,1:NJ),Ap(i,1:NJ),w*An(i,1:NJ),RHSj(1:NJ));
        for i=2:NI-1
%             RHSj(2:NJ-1) = (1-w)*Ap(i,2:NJ-1).*phi(i,2:NJ-1,2)-w*Ae(i,2:NJ-1).*phi(i+1,2:NJ-1,2)-w*Aw(i,2:NJ-1).*phi(i-1,2:NJ-1,3)+w*S(i,2:NJ-1,t);
%             RHSj(1 )     = (1-w)*Ap(i,     1).*phi(i,     1,2)-w*Ae(i,     1).*phi(i+1,     1,2)-w*Aw(i,     1).*phi(i-1,     1,3)+w*g(i,     1,t);   %south BC
%             RHSj(NJ)     = (1-w)*Ap(i,    NJ).*phi(i,    NJ,2)-w*Ae(i,    NJ).*phi(i+1,    NJ,2)-w*Aw(i,    NJ).*phi(i-1,    NJ,3)+w*g(i,    NJ,t);   %north BC
            RHSj(1:NJ) = (1-w)*Ap(i,1:NJ).*phi(i,1:NJ,2)-w*Ae(i,1:NJ).*phi(i+1,1:NJ,2)-w*Aw(i,1:NJ).*phi(i-1,1:NJ,3)+w*S(i,1:NJ,t);
            %solve for k+1
            phi(i,1:NJ,3) = TDS1(1,NJ,w*As(i,1:NJ),Ap(i,1:NJ),w*An(i,1:NJ),RHSj(1:NJ));
        end
        i=NI;    %east
%          RHSj(2:NJ-1) = (1-w)*Ap(i,2:NJ-1).*phi(i,2:NJ-1,2)-w*Aw(i,2:NJ-1).*phi(i-1,2:NJ-1,3)+w*f(i,2:NJ-1,t);
%          RHSj(1 )     = (1-w)*Ap(i,     1).*phi(i,     1,2)-w*Aw(i,     1).*phi(i-1,     1,3)+w*g(i,     1,t);   %south BC
%          RHSj(NJ)     = (1-w)*Ap(i,    NJ).*phi(i,    NJ,2)-w*Aw(i,    NJ).*phi(i-1,    NJ,3)+w*f(i,    NJ,t);   %north BC
         RHSj(1:NJ) = (1-w)*Ap(i,1:NJ).*phi(i,1:NJ,2)-w*Aw(i,1:NJ).*phi(i-1,1:NJ,3)+w*S(i,1:NJ,t);
         %solve for k+1
         phi(i,1:NJ,3) = TDS1(1,NJ,w*As(i,1:NJ),Ap(i,1:NJ),w*An(i,1:NJ),RHSj(1:NJ));
        
        %find the rms residual and errors
        rho_tot = 0;
        e_tot = 0;
        for j=1:NJ
            for i=1:NI
                rho_tot = rho_tot + (Ap(i,j)*phi(i,j,3)-Ap(i,j)*phi(i,j,1))^2; 
                e_tot   =  e_tot  + (        phi(i,j,3)-        phi(i,j,1))^2;
            end
        end
        rho_rms(k) = sqrt(rho_tot/((NI)*(NJ)));
        e_rms(k)   = sqrt( e_tot /((NI)*(NJ)));
        
        resid  = rho_rms(k);
        lambda = e_rms(k)/e_rms(k-1);
        
        if rem(k,loop_iterations)==0
            if lambda > 1.01 & w > 0.01;
                w = w - 0.01;
            elseif lambda > 0.99995 & w < 0.1
                w = w + 0.005;
            end
            fprintf('  %4i  %6i  %11.4g  %11.4g  %5.3f  %10.3f\n',t,k,resid,lambda,w,toc)
            %disp([t k w resid lambda])
            figure(941),semilogy(rho_rms);
%             pause
        %plot current result for testing
        figure(944),surf(x,y,phi(:,:,3)),pause(0.0001)
        end
        

%         %plot current result for testing
%         figure(944),surf(x,y,phi(:,:,3)),pause(0.0001)
        
    end %solution is converged
    
    %store final solution into pressure table
    p(:,:,t) = phi(:,:,3);
    
    %convergence information
%     t
%     k
%     resid
%     lambda
%     w
%    toc
    fprintf('  %4i  %6i  %11.4g  %11.4g  %5.3f  %10.3f\n',t,k,resid,lambda,w,toc)
%     if k>1
%         lamda_max=e_rms(k)/e_rms(k-1)
%     end
    
    %plot convergence
    figure(942),hold on,semilogy(ktot+[2:k],rho_rms(2:k)),title(['rho_r_m_s for timestep t=',num2str(t)])
    hold off
    %figure(943),plot(rho_rms),title(['e_r_m_s for timestep t=',num2str(t)])
    %pause
       
    ktot = ktot + k;
    
%     %use an offset to peg the reference pressure over time... 
%     %is this okay to do?
%     if nargin>=9
%         %pause
%         deltaP = Pref(t) - p(Px,Py,t);
%         p(:,:,t) = p(:,:,t) + deltaP * ones(NI,NJ);
%     end

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
