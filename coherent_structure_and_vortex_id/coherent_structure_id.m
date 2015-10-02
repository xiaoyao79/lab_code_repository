function [cal_mat] = coherent_structure_id(u,v,dx,dy,mask,vortex_method)
[NY,NX,NT] = size(u);

if size(mask,3) ~= size(u,3)
    MASK = repmat(mask,[1 1 size(u,3)]);
else
    MASK = mask;
end

[dudx] = socdiff_bc(u,dx,2,MASK);
[dudy] = socdiff_bc(u,dx,1,MASK);
[dvdx] = socdiff_bc(v,dy,2,MASK);
[dvdy] = socdiff_bc(v,dy,1,MASK);
vort   = dvdx - dudy;

cal_mat = zeros([NY NX NT]);
swirl   = zeros([NY NX NT]);
stretch = zeros([NY NX NT]);

for t=1:NT
    for i=1:NX
        for j=1:NY
            vgrad = [dudx(j,i,t),dudy(j,i,t);dvdx(j,i,t),dvdy(j,i,t)];
            Omega = 0.5*(vgrad+transpose(vgrad));
            S     = 0.5*(vgrad-transpose(vgrad));  
            if strcmp(vortex_method,'Qcrit')
                Omegan=norm(Omega,'fro');
                Sn=norm(S,'fro');
                Q=0.5*(Omegan.^2-Sn.^2);
                if Q>0
                    cal_mat(j,i,t)=Q;
                end
            elseif strcmp(vortex_method,'Dcrit')
                Omegan=norm(Omega,'fro');
                Sn=norm(S,'fro');
                Q=0.5*(Omegan.^2-Sn.^2);
                Delta=(Q/3)^3+(det(vgrad)/2)^2;
                if Delta>0
                    cal_mat(j,i,t)=Delta;
                end
            elseif strcmp(vortex_method,'Lambda2')
                T=eig(S^2+Omega^2);
%                 Lamda_min = min(T);
                Lamda_max = max(T);
                if Lamda_max < 0;
                    cal_mat(j,i,t)=Lamda_max;
                else
                    cal_mat(j,i,t)=0;
                end
            elseif strcmp(vortex_method,'LambdaCI')
                EV = eig(vgrad);
                if ~isreal(EV(1))
                    swirl(j,i,t)   = abs(imag(EV(1)));
                    stretch(j,i,t) = real(EV(1));
                end
            end
        end
    end
end

if strcmp(vortex_method,'LambdaCI')
    cal_mat = swirl./max(swirl(:)).*sign(vort);
else
    cal_mat = cal_mat.*sign(vort);
end

