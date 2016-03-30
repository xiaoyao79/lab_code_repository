function [Upsi,Vpsi,psi,omegapsi] = cdev_cartesian(U,V,dx,dy,omega,mask)
% PSI-OMEGA SOLVER Iterative solver for vector field reconstruction
ave_weight_case = 'weighted';  % Either 'unweighted' or 'weighted'
solver_case     = 'velocity'; % Either 'vorticity' or 'velocity'

% Determine physical domain size
[s1,s2] = size(U);

% Discrete Psi
% Compute boundary Psi values (boundary conditions for solver)
% Initialize Psi Matrix for boundary value storage
psi      = zeros(s1,s2);

% Calculate boundary values across top of domain
psi(1,:) = -cumtrapz(V(1,:),2)*dx;

% Calculate boundary values along right side of domain
if max(abs(U(:))) ~= 0
    psi(:,end)  = cumtrapz(U(:,end),1)*dy + psi(1,end);
else
    intstart = psi(:,end);
    intdiff = 1;
    while intdiff >1E-7
        for k = 2:size(V,1)-1
            psi(k,end) = 0.5*(psi(k-1,end)+psi(k+1,end))+0.25*omega(k,end)*dy^2;
        end
        intdiff = sum((1/numel(intstart))*sqrt((intstart-psi(:,end)).^2));
        intstart = psi(:,end);
    end
end

% Calculate boundary values across bottom of domain
temp0            = fliplr(cumtrapz(fliplr(V(end,:)),2))*dx+psi(end,end);
psi(end,1:end-1) = temp0(1:end-1);

% Calculate boundary values along left side of domain
if max(abs(U(:))) ~= 0
    psi(:,1) = -flipud(cumtrapz(flipud(U(:,1)),1))*dy + psi(end,1);
else
    intstart = psi(:,1);
    intdiff = 1;
    while intdiff >1E-7
        for k = 2:size(V,1)-1
            psi(k,1) = 0.5*(psi(k-1,1)+psi(k+1,1))+0.25*omega(k,1)*dy^2;
        end
        intdiff = sum((1/numel(intstart))*sqrt((intstart-psi(:,end)).^2));
        intstart = psi(:,1);
    end
end

% Integrate Psi values from Left-to-Right across domain
psil = zeros(s1,s2-1);
for i = 1:s2-1
    if i == 1
        psil(:,i) = -0.5*(V(:,2)+V(:,1))*dx;
    elseif i == 2
        psil(:,i) = -0.5*(V(:,3)+2*V(:,2)+V(:,1))*dx;
    else
        psil(:,i) = -0.5*(V(:,i+1)+2*sum(V(:,2:i),2)+V(:,1))*dx;
    end
end
psil = psil + repmat(psi(:,1),[1 size(psil,2)]);
psil = cat(2,psi(:,1),psil);

% Integrate Psi values from Right-to-Left across domain
psir  = zeros(s1,s2-1);
Vflip = fliplr(-V);
for i = 1:s2-1
    if i == 1
        psir(:,i) = -0.5*(Vflip(:,2)+Vflip(:,1))*dx;
    elseif i == 2
        psir(:,i) = -0.5*(Vflip(:,3)+2*Vflip(:,2)+Vflip(:,1))*dx;
    else
        psir(:,i) = -0.5*(Vflip(:,i+1)+2*sum(Vflip(:,2:i),2)+Vflip(:,1))*dx;
    end
end
psir = psir + repmat(psi(:,end),[1 size(psir,2)]);
psir = (cat(2,psi(:,end),psir));

% Integrate Psi values from Left-to-Right across domain
psiu = zeros(s1-1,s2);
for i = 1:s1-1
    if i == 1
        psiu(i,:) = -0.5*(U(2,:)+U(1,:))*dx;
    elseif i == 2
        psiu(i,:) = -0.5*(U(3,:)+2*U(2,:)+U(1,:))*dx;
    else
        psiu(i,:) = -0.5*(U(i+1,:)+2*sum(U(2:i,:),1)+U(1,:))*dx;
    end
end
psiu = psiu + repmat(psi(1,:),[size(psiu,1) 1]);
psiu = cat(1,psi(1,:),psiu);

% Integrate Psi values from Right-to-Left across domain
psid  = zeros(s1-1,s2);
Uflip = flipud(-U);
for i = 1:s1-1
    if i == 1
        psid(i,:) = -0.5*(Uflip(2,:)+Uflip(1,:))*dx;
    elseif i == 2
        psid(i,:) = -0.5*(Uflip(3,:)+2*Uflip(2,:)+Uflip(1,:))*dx;
    else
        psid(i,:) = -0.5*(Uflip(i+1,:)+2*sum(Uflip(2:i,:),1)+Uflip(1,:))*dx;
    end
end
psid = psid + repmat(psi(end,:),[size(psid,1) 1]);
psid = (cat(1,psi(end,:),psid));

% Stream Function Computation
switch ave_weight_case
    case 'unweighted'
        % Calculate Psi based on averaging
        psi = (psil+fliplr(psir))/2;
    case 'weighted'
        % Calculate Psi based on weighted averaging
        weight_right  = repmat(linspace(0,1,size(psi,2)),[size(psi,1) 1]);
        weight_left   = repmat(linspace(1,0,size(psi,2)),[size(psi,1) 1]);
        
        psi = (weight_left.*psil + weight_right.*fliplr(psir));
    otherwise
        warning('Unexpected averaging setting.');
end

% Set initial count value (n), Error value (SSE), and initial Psi estimate
SSE  = 1;
n    = 0;
psi1 = psi;

switch solver_case
    case 'vorticity'
        while SSE >= 1E-5
            % Compute Updated Psi from Vorticity & Psi Field
            for i = 2:s1-1
                for j = 2:s2-1
                    if mask(i,j) == 1
                        psi(i,j) = 0.25*(psi(i+1,j)... % Right of the current cell
                            +psi(i-1,j)...             % Left of the current cell
                            +psi(i,j+1)...             % Above current cell
                            +psi(i,j-1)...             % Below current cell
                            +omega(i,j)*dx^2);         % vorticity in current cell
                    end
                end
            end
            
            omegapsi = -(socdiff2(psi,dx,2)+socdiff2(psi,dx,1));
            
            SSE      = omega-omegapsi;
            SSE      = sqrt(mean(SSE(:).^2))/(max(omega(:))-min(omega(:)));
            
            Vpsi = -socdiff(psi,dx,2);
            Upsi =  socdiff(psi,dx,1);
            
            n = n+1;
            
            TREND(n) = SSE;
            
            figure(100); plot(TREND);
            set(gca,'YScale','log');
            
            if n > 2
                if abs(TREND(n)-TREND(n-1)) < 1E-20
                    break
                end
                if n == 100
                    break
                end
            end
        end
    case 'velocity'
        while SSE >= 1E-5
            % Compute Updated Psi from Vorticity & Psi Field
            for i = 2:s1-1
                for j = 2:s2-1
                    if mask(i,j) == 1
                        psi(i,j) = 0.25*(psi(i+1,j)... % Right of the current cell
                            +psi(i-1,j)...             % Left of the current cell
                            +psi(i,j+1)...             % Above current cell
                            +psi(i,j-1)...             % Below current cell
                            +omega(i,j)*dx^2);         % vorticity in current cell
                    end
                end
            end
            
            Vpsi =  V;
            Upsi =  socdiff_bc(psi,dx,1,mask);
            
            SSE = norm(psi1-psi,Inf)/norm(psi1,Inf);
            
            n = n+1;
            
            TREND(n) = SSE;
            
            figure(100); plot(TREND);
            set(gca,'YScale','log');
            pause(5E-2)
            
            if n > 2
                if abs(TREND(n)-TREND(n-1)) < 1E-20
                    break
                end
                if n == 2000
                    break
                end
            end
            
            omega = socdiff_bc(Vpsi,dx,2,mask) - socdiff_bc(Upsi,dx,1,mask);
            omega = laplace_outlier_detect(omega,dx,dy,1);
            
            psi1 = psi;
        end
end
omegapsi = omega;
end