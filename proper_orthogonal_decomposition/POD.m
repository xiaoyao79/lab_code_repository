function [Us,Vs] = POD(U,V)
pc_NX = size(U,1);
pc_NY = size(U,2);

Umean = mean(U,3);
Vmean = mean(V,3);
Us = U;
Vs = V;

for t=1:size(U,3)
    Us(:,:,t) = U(:,:,t)-Umean;
    Vs(:,:,t) = V(:,:,t)-Vmean;
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

for t=1:min(size(U,3),40)
    figure(18)
    subplot(10,4,t)
    plot(A(t,:))
    ylabel(num2str(t))
end


%rebuild POD modes
NN = find(cumsum(D)/sum(D)>.95,1,'first');
NN = max(NN,10);

%moving average smoothing
As = A;
As(:,1        ) = 1/5*(A(:,1        ) + A(:,2        ) + A(:,3        ) + A(:,  size(U,3)-1) + A(:,  size(U,3)-0) );
As(:,2        ) = 1/5*(A(:,1        ) + A(:,2        ) + A(:,3        ) + A(:,4        ) + A(:,  size(U,3)-0) );
As(:,3:size(U,3)-2) = 1/5*(A(:,1:size(U,3)-4) + A(:,2:size(U,3)-3) + A(:,3:size(U,3)-2) + A(:,4:size(U,3)-1) + A(:,5:size(U,3)-0) );
As(:,  size(U,3)-1) = 1/5*(A(:,1        ) + A(:,  size(U,3)-3) + A(:,  size(U,3)-2) + A(:,  size(U,3)-1) + A(:,  size(U,3)-0) );
As(:,  size(U,3)  ) = 1/5*(A(:,1        ) + A(:,2        ) + A(:,  size(U,3)-2) + A(:,  size(U,3)-1) + A(:,  size(U,3)-0) );

for t=1:min(size(U,3),40)
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
% keyboard
Us = repmat(Umean,[1,1,size(U,3)])+reshape(reshape(M1(:,:,1:NN),pc_NX*pc_NY,NN)*As(1:NN,:),pc_NX,pc_NY,size(U,3));
Vs = repmat(Vmean,[1,1,size(U,3)])+reshape(reshape(M2(:,:,1:NN),pc_NX*pc_NY,NN)*As(1:NN,:),pc_NX,pc_NY,size(U,3));

fprintf('%g\n',toc)

Us = Us;
Vs = Vs;