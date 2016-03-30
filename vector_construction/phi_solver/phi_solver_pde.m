function [Uphi,Vphi,phi,Qphi] = phi_solver_pde(U,V,dx,dy,Q,mask)
% PSI-OMEGA SOLVER Iterative solver for vector field reconstruction
ave_weight_case = 'weighted';  % Either 'unweighted' or 'weighted'

% Determine physical domain size
[s1,s2] = size(U);

% Discrete Psi
% Compute boundary Psi values (boundary conditions for solver)
% Initialize Psi Matrix for boundary value storage
phi = zeros(s1,s2);

intstart = vertcat(phi(1:end,1),phi(end,2:end-1)',phi(2:end-1,end),phi(1,2:end)');
intdiff = 1;
while intdiff > 1E-7
    for k = 2:size(Q,2)-1
        phi(1,k) = 0.5*(phi(1,k-1)+phi(1,k+1)+Q(1,k)*dx^2);
    end
    phi(1,end) = 0.5*(phi(1,end-1)+phi(2,end)+Q(1,end)*dx^2);
    for k = 2:size(Q,1)-1
        phi(k,end) = 0.5*(phi(k-1,end)+phi(k+1,end)+Q(k,end)*dx^2);
    end
    phi(end,end) = 0.5*(phi(end-1,end)+phi(end,end-1)+Q(end,end)*dx^2);
    for k = 2:size(Q,1)-1
        phi(end,k) = 0.5*(phi(end,k-1)+phi(end,k+1)+Q(end,k)*dx^2);
    end
    phi(end,1) = 0.5*(phi(end,2)+phi(end-1,1)+Q(end,1)*dx^2);
    for k = 2:size(Q,1)-1
        phi(k,1) = 0.5*(phi(k-1,1)+phi(k+1,1)+Q(k,1)*dx^2);
    end
    intcur  = vertcat(phi(1:end,1),phi(end,2:end-1)',phi(2:end-1,end),phi(1,2:end)');
    intdiff = sum((1/numel(intstart))*sqrt((intstart-intcur).^2));
    intstart = intcur;
end

SSE = 1;
n   = 0;

phi1 = phi;

while SSE >= 1E-5
    % Compute Updated Psi from Vorticity & Psi Field
    for i = 2:s1-1
        for j = 2:s2-1
            if mask(i,j) == 1
                phi(i,j) = 0.25*(phi(i+1,j)... % Right of the current cell
                    +phi(i-1,j)...             % Left of the current cell
                    +phi(i,j+1)...             % Above current cell
                    +phi(i,j-1)...             % Below current cell
                    +Q(i,j)*dx^2);         % vorticity in current cell
            end
        end
    end
    
    Vphi =  socdiff_bc(phi,dx,1,mask);
    Uphi =  socdiff_bc(phi,dx,2,mask);
    Qphi = Q;

    SSE = norm(phi1-phi,Inf)/norm(phi1,Inf);
    
    n = n+1;
    TREND(n) = SSE;
    
    figure(100); plot(TREND);
    set(gca,'YScale','log');
    pause(5E-2);
    
    phi1 = phi;
end
end