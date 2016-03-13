function [Us,Vs,mask] = run_POD(Ut,Vt,mask,pc_NX,pc_NY,penergy,resamp,smoothtrig)
%%
Umean = mean(Ut.*mask,3);
Vmean = mean(Vt.*mask,3);

Us = Ut;
Vs = Vt;

for t=1:size(Vt,3)
    Us(:,:,t) = Ut(:,:,t).*mask(:,:,t)-Umean;
    Vs(:,:,t) = Vt(:,:,t).*mask(:,:,t)-Vmean;
end

tic,fprintf('POD filtering velocity errors...')

[M1,M2,D,A] = KLDs2v(Us,Vs);

figure(17)
subplot(2,1,1)
loglog(D./sum(D),'-o')
title('Eigenspectra')
ylabel('Eigenvalue'); xlabel('Mode')
ylim([10^-10 10^0])
subplot(2,1,2)
semilogx(cumsum(D)./sum(D),'-o')
title('Psuedo Energy')
ylabel('Percent Energy'); xlabel('Mode')

for t=1:min(size(Vs,3),40)
    figure(18)
    subplot(10,4,t)
    plot(A(t,:))
    ylabel(num2str(t))
end

%rebuild POD modes
NN = find(cumsum(D)/sum(D)>penergy,1,'first');
NN = max(NN,10);
size(Vs,3)

% Super-sample Coefficient Matrix
[xgrid0,ygrid0] = meshgrid(linspace(1,size(A,1),size(A,1)),...
    linspace(1,size(A,2),size(A,2)));
[xgrid1,ygrid1] = meshgrid(linspace(1,size(A,1),resamp*size(A,1)),...
    linspace(1,size(A,2),resamp*size(A,2)));

A = interp2(xgrid0,ygrid0,A,xgrid1,ygrid1,'cubic');

%moving average smoothing
As = A;
if strcmp(smoothtrig,'on')
    As(:,1             ) = 1/5*(A(:,1        )      + A(:,2        )      + A(:,3        )      + A(:,  size(Vs,3)-1) + A(:,  size(Vs,3)-0) );
    As(:,2             ) = 1/5*(A(:,1        )      + A(:,2        )      + A(:,3        )      + A(:,4        )      + A(:,  size(Vs,3)-0) );
    As(:,3:size(Vs,3)-2) = 1/5*(A(:,1:size(Vs,3)-4) + A(:,2:size(Vs,3)-3) + A(:,3:size(Vs,3)-2) + A(:,4:size(Vs,3)-1) + A(:,5:size(Vs,3)-0) );
    As(:,  size(Vs,3)-1) = 1/5*(A(:,1        )      + A(:,  size(Vs,3)-3) + A(:,  size(Vs,3)-2) + A(:,  size(Vs,3)-1) + A(:,  size(Vs,3)-0) );
    As(:,  size(Vs,3)  ) = 1/5*(A(:,1        )      + A(:,2        )      + A(:,  size(Vs,3)-2) + A(:,  size(Vs,3)-1) + A(:,  size(Vs,3)-0) );
end

for t=1:min(resamp*size(Vs,3),40)
    figure(18)
    subplot(10,4,t)
    hold on
    plot(As(t,:),'r')
    ylabel(num2str(t))
    hold off
end

fprintf('reconstructing %u modes',NN)
%NN = 1:NN;
fprintf(' and %.6g percent energy\n', sum(D(1:NN))/sum(D)*100 )

%use NN timesteps with smoothed coefficients
Us = repmat(Umean,[1,1,resamp*size(Vs,3)]) + reshape(reshape(M1(:,:,1:NN), pc_NX*pc_NY,NN) * As(1:NN,:), pc_NX,pc_NY,resamp*size(Vs,3));
Vs = repmat(Vmean,[1,1,resamp*size(Vs,3)]) + reshape(reshape(M2(:,:,1:NN), pc_NX*pc_NY,NN) * As(1:NN,:), pc_NX,pc_NY,resamp*size(Vs,3));

fprintf('%g\n',toc)

h = fspecial('average',[5 5]);

[xgrid2,ygrid2,zgrid2] = meshgrid(linspace(1,size(Us,2),size(Us,2)),...
    linspace(1,size(Us,1),size(Us,1)),linspace(1,size(Ut,3),size(Ut,3)));
[xgrid3,ygrid3,zgrid3] = meshgrid(linspace(1,size(Us,2),size(Us,2)),...
    linspace(1,size(Us,1),size(Us,1)),linspace(1,size(Ut,3),size(Us,3)));

% keyboard

mask = interp3(xgrid2,ygrid2,zgrid2,mask,xgrid3,ygrid3,zgrid3);
mask(mask~=0) = 1;

for t = 1:size(Us,3)
    Us(:,:,t) = imfilter(Us(:,:,t),h).*mask(:,:,t);
    Vs(:,:,t) = imfilter(Vs(:,:,t),h).*mask(:,:,t);
end