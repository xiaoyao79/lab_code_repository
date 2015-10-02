function [Us,Vs] = run_POD_ave_filt(Useries1,Useries2,Vseries1,Vseries2,mask,pc_NX,pc_NY,rep_switch)
%  keyboard
Umean1 = mean(Useries1.*mask,3);
Vmean1 = mean(Vseries1.*mask,3);

Umean2 = mean(Useries2.*mask,3);
Vmean2 = mean(Vseries2.*mask,3);

Us1 = Useries1;
Vs1 = Vseries1;

Us2 = Useries2;
Vs2 = Vseries2;

if rep_switch == 1
    Umean = mean([Umean1(:) Umean2(:)]');
    Umean = reshape(Umean,[size(Umean1,1) size(Umean1,2) size(Umean1,3)]);
    
    Vmean = mean([Vmean1(:) Vmean2(:)]');
    Vmean = reshape(Vmean,[size(Vmean1,1) size(Vmean1,2) size(Vmean1,3)]);
else
    Umean = Umean1;
    Vmean = Vmean1;
end

%%
for t=1:size(Vseries1,3)
    Us1(:,:,t) = Useries1(:,:,t).*mask(:,:,t)-Umean1;
    Vs1(:,:,t) = Vseries1(:,:,t).*mask(:,:,t)-Vmean1;

    Us2(:,:,t) = Useries2(:,:,t).*mask(:,:,t)-Umean2;
    Vs2(:,:,t) = Vseries2(:,:,t).*mask(:,:,t)-Vmean2;
end
%%
tic,fprintf('POD filtering velocity errors...')

[M11,M21,D1,A1] = KLDs2v(Us1,Vs1);
[M12,M22,D2,A2] = KLDs2v(Us2,Vs2);
%%
if rep_switch == 1
for t = 1:size(M11,3)
    M1(:,:,t) = (M11(:,:,t) + M12(:,:,t)) / 2;
    M2(:,:,t) = (M21(:,:,t) + M22(:,:,t)) / 2;
end
else
    M1 = M11;
    M2 = M21;
end

D  = (D1 + D2) / 2;
A  = (A1 + A2) / 2;

%%
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


for t=1:min(size(Vs1,3),40)
    figure(18)
    subplot(10,4,t)
    plot(A(t,:))
    ylabel(num2str(t))
end

%rebuild POD modes
NN = find(cumsum(D)/sum(D)>.90,1,'first');
NN = max(NN,10);
size(Vs1,3)

%moving average smoothing
As = A;
As(:,1              ) = 1/5*(A(:,1        )       + A(:,2        )       + A(:,3        )       + A(:,  size(Vs1,3)-1) + A(:,  size(Vs1,3)-0) );
As(:,2              ) = 1/5*(A(:,1        )       + A(:,2        )       + A(:,3        )       + A(:,4        )       + A(:,  size(Vs1,3)-0) );
As(:,3:size(Vs1,3)-2) = 1/5*(A(:,1:size(Vs1,3)-4) + A(:,2:size(Vs1,3)-3) + A(:,3:size(Vs1,3)-2) + A(:,4:size(Vs1,3)-1) + A(:,5:size(Vs1,3)-0) );
As(:,  size(Vs1,3)-1) = 1/5*(A(:,1        )       + A(:,  size(Vs1,3)-3) + A(:,  size(Vs1,3)-2) + A(:,  size(Vs1,3)-1) + A(:,  size(Vs1,3)-0) );
As(:,  size(Vs1,3)  ) = 1/5*(A(:,1        )       + A(:,2        )       + A(:,  size(Vs1,3)-2) + A(:,  size(Vs1,3)-1) + A(:,  size(Vs1,3)-0) );

for t=1:min(size(Vs1,3),40)
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
Us = repmat(Umean,[1,1,size(Vs1,3)]) + reshape(reshape(M1(:,:,1:NN), pc_NX*pc_NY,NN) * As(1:NN,:), pc_NX,pc_NY,size(Vs1,3));
Vs = repmat(Vmean,[1,1,size(Vs1,3)]) + reshape(reshape(M2(:,:,1:NN), pc_NX*pc_NY,NN) * As(1:NN,:), pc_NX,pc_NY,size(Vs1,3));
%%

fprintf('%g\n',toc)
% keyboard
Us = Us.*mask;
Vs = Vs.*mask;
