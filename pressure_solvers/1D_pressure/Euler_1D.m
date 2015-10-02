function [dPdx,P,ux,ut] = Euler_1D(U,DX,DT)
% Euler_1D Pressure determination from 1D+t velocity fields
%
% [dPdx,P,u_x,u_t] = Euler_1D(U,DX,DT,MU,RHO) requires inputs for the 1D+t
% velocity field, U (in matrix form), and physical and temporal spacing,
% DX & DT.

NU  = 3.77e-6;  % Kinematic Viscosity, m^2/s
RHO = 1100;     % Density            , kg/m^3
MU  = NU*RHO;   % Dynamic Viscosity  , kg/m-s

dudt   = socdiff (U,DT,2);
dudx   = socdiff (U,DX,1);
d2udx2 = socdiff2(U,DX,1);

dPdx  = MU*(d2udx2) - RHO*(dudt+U.*dudx);

P = zeros(size(dPdx));

for i=3:size(dPdx,1)
    P(i,:) = simpsonH(dPdx(1:i,:),DX);
end

ux = U.^2/2;
ut = zeros(size(dudt));

for i=3:size(dudt,1)
    ut(i,:) = simpsonH(dudt(1:i,:),DX);
end