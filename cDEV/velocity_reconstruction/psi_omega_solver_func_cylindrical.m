function [Utpsi,Urpsi,psi,omegapsi] = ...
    psi_omega_solver_func_cylindrical(Ut,Ur,r,dr,dt,omega,mask,num_iter)
% PSI-OMEGA SOLVER Iterative solver for vector field reconstruction

ave_weight_case = 'weighted';  % Either 'unweighted' or 'weighted'
solver_case     = 'velocity';  % Either 'vorticity' or 'velocity'

% Determine physical domain size
[s1,s2] = size(Ut);

% Discrete Psi
% Compute boundary Psi values (boundary conditions for solver)

% Initialize Psi Matrix for boundary value storage
psi = zeros(s1,s2);

% Calculate boundary values across top of domain
psi(1,:) = r(1,:).*cumtrapz(Ur(1,:),2)*dt;
psit     = r(1,:).*(cumtrapz(fliplr(Ur(1,:)),2))*dt;

% Calculate boundary values along right side of domain
if max(Ut(:)) ~= 0
    temp0  = -cumtrapz(Ut(:,end),1)*dr + psi(1,end);
else
    temp0  = -cumtrapz(Ut(:,end),1)*dr + psi(1,end);
    %     temp0  = cumtrapz(cumtrapz(omega(:,end)/4,1),1)*dr^2 + psi(1,end);
end
psi(:,end) = temp0;

% Calculate boundary values across bottom of domain
temp0            = -r(end,:).*fliplr(cumtrapz(fliplr(Ur(end,:)),2))*dt + psi(end,end);
% temp0            = -r(end,:).*fliplr(cumtrapz(fliplr(Ur(end,:)),2))*dt;
psi(end,1:end-1) = temp0(1:end-1);

% Calculate boundary values along left side of domain
if max(Ut(:)) ~= 0
    temp0 = flipud(cumtrapz(flipud(Ut(:,1)),1))*dr + psi(end,1);
    %     temp0 = -flipud(cumtrapz(flipud(Ut(:,1)),1))*dr + psit(1,1);
else
    temp0 = flipud(cumtrapz(flipud(Ut(:,1)),1))*dr + psi(end,1);
    %     temp0 = -flipud(cumtrapz(cumtrapz(flipud(-omega(:,1)/4),1),1))*dr^2 +...
    %         psi(end,1);
end
psi(:,1)  = temp0;

% Integrate Psi values from Left-to-Right across domain
psil = zeros(s1,s2-1);

for i = 1:s2-1
    if i == 1
        psil(:,i) = 0.5*r(:,i).*(Ur(:,2)+Ur(:,1))*dt;
    elseif i == 2
        psil(:,i) = 0.5*r(:,i).*(Ur(:,3)+2*Ur(:,2)+Ur(:,1))*dt;
    else
        psil(:,i) = 0.5*r(:,i).*(Ur(:,i+1)+2*sum(Ur(:,2:i),2)+Ur(:,1))*dt;
    end
end

psil = psil + repmat(psi(:,1),[1 size(psil,2)]);
psil = cat(2,psi(:,1),psil);

% Integrate Psi values from Right-to-Left across domain
psir  = zeros(s1,s2-1);

Vflip = fliplr(-Ur);

for i = 1:s2-1
    if i == 1
        psir(:,i) = 0.5*r(:,i).*(Vflip(:,2)+Vflip(:,1))*dt;
    elseif i == 2
        psir(:,i) = 0.5*r(:,i).*(Vflip(:,3)+2*Vflip(:,2)+Vflip(:,1))*dt;
    else
        psir(:,i) = 0.5*r(:,i).*(Vflip(:,i+1)+2*sum(Vflip(:,2:i),2)+Vflip(:,1))*dt;
    end
end

psir = psir + repmat(psi(:,end),[1 size(psir,2)]);
psir = (cat(2,psi(:,end),psir));

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
        % psi(2:end-1,2:end-1) = temp1(2:end-1,2:end-1);
    otherwise
        warning('Unexpected averaging setting.');
end

SSE = 1;
n   = 0;

psi1 = psi;

Urpsi = Ur;
Utpsi = -socdiff(psi,dr,1);

% figure(1000);
% subplot(2,2,1);
% imagesc(Urpsi); colorbar;
% subplot(2,2,2);
% imagesc(Utpsi); colorbar;
% subplot(2,2,3);
% imagesc(psi); colorbar;
% subplot(2,2,4);
% imagesc(omega); colorbar;

% Generate boundary mask for use in central differencing codes
bmask = mask;
bmask(bmask ~= 0) = 1;

% Paused here 22:27 on 07/25/2015
switch solver_case
    case 'vorticity'
        while SSE >= 1E-4
            % Compute Updated Psi from Vorticity & Psi Field
            for i = 2:s1-1
                for j = 2:s2-1
                    scale = (2*r.^2*dt^2+2*dr^2)^-1;
                    psi(i,j) = scale.*(r(i,j)*(psi(i+1,j)-psi(i-1,j))*dr*dt^2 ...
                        +r(i,j)^2*(psi(i+1,j)+psi(i-1,j))*dt^2 ...
                        +(psi(i,j+1)+psi(i,j-1)) ...
                        +r(i,j)^2*omega(i,j)*dr^2*dt^2);
                end
            end
            
            omegapsi = (r.^-1).*(Utpsi+ socdiff(Utpsi,dr,1) - socdiff(Urpsi,dt,2));
            
            SSE      = omega-omegapsi;
            SSE      = sqrt(mean(SSE(:).^2))/(max(omega(:))-min(omega(:)));
            
            Urpsi = (r.^-1).*socdiff(psi,dt,2);
            Utpsi =  -socdiff(psi,dr,1);
            
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
        while SSE >= eps
            % Compute Updated Psi from Vorticity & Psi Field
            for i = 2:s1-1
                for j = 2:s2-1
                    if mask(i,j) > 1
                        scale = (4*r.^2*dt^2+4*dr^2).^-1;
                        psi(i,j) = scale(i,j)*(r(i,j)*(psi(i+1,j)-psi(i-1,j))*dr*dt^2 ...
                            +              2*r(i,j)^2*(psi(i+1,j)+psi(i-1,j))*dt^2 ...
                            +                       2*(psi(i,j+1)+psi(i,j-1))*dr^2 ...
                            +2*r(i,j)^2*omega(i,j)*dr^2*dt^2);
                    end
                end
            end
            
            if min(mask(:)) == 1
                Urpsi = (r.^-1).*socdiff(psi,dt,2);
            else
                Urpsi = (r.^-1).*socdiff_bc(psi,dt,2,bmask);
            end
            q = [1 1];
            
            %             figure(1000);
            %             subplot(2,2,1);
            %             imagesc(Urpsi); colorbar;
            
            for i = 1:s1
                for j = 1:s2
                    if abs(Ur(i,j)-Urpsi(i,j)) > 0.05*abs(Ur(i,j)) && mask(i,j) > 2
                        Urpsi(i,j) = Ur(i,j);
                    end
                end
            end
            
%                         Ur = Urpsi;
%                         Urpsi = Ur;

            if n == 0
                Ut_old = Ut;
            else
                Ut_old = Utpsi;
            end
            
            if min(mask(:)) == 1
                Utpsi =  -socdiff(psi,dr,1);
            else
                Utpsi =  -socdiff_bc(psi,dr,1,bmask);
            end
            
            if n == 0
                omegapsi_old = omega;
            else
                omegapsi_old = omegapsi;
            end
            
            if min(mask(:)) == 1
                omegapsi = (r.^-1).*...
                    (socdiff(r.*Utpsi,dr,1) - socdiff(Urpsi,dt,2));
            else
                omegapsi = (r.^-1).*...
                    (socdiff_bc(r.*Utpsi,dr,1,bmask) -...
                    socdiff_bc(Urpsi,dt,2,bmask));
            end
            
            omegapsi = laplace_outlier_detect(omegapsi,dt,dr,2);
            
            fprintf('   Pass %3.0f Residual %04E \n',n,SSE);
            
%             SSE = norm(omegapsi-omegapsi_old)/norm(omegapsi_old);
%             SSE = norm(Utpsi-Ut_old,'inf');
            SSE = rms(Utpsi(:)-Ut_old(:));
%             keyboard
            n = n+1;
            
            TREND(n) = SSE;
            
            %             figure(100); plot(TREND);
            %             set(gca,'YScale','log');
            
            if n > 2
%                 if abs(TREND(n)-TREND(n-1)) < eps
                if TREND(n) < eps
                    break
                end
                if n == num_iter
                    break
                end
            end
            
                        figure(1000);
                        subplot(2,2,1);
                        imagesc(Urpsi); colorbar;
                        subplot(2,2,2);
                        imagesc(Utpsi); colorbar;
                        subplot(2,2,3);
                        imagesc(psi); colorbar;
                        subplot(2,2,4);
                        imagesc(omegapsi); colorbar;
            %
            %             omega = (r.^-1).*(socdiff(r.*Utpsi,dr,1) - socdiff(Urpsi,dt,2));
            %             omega = laplace_outlier_detect(omega,dt,dr,1);
            omega = omegapsi;
%             keyboard
        end
end
omegapsi = omega;
end