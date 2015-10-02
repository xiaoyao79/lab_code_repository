function [Displx_new,Disply_new] = tss(Displx,Disply,facq,SSF,grid_scale)

% The routine input is given by two displacement vector fields
% (two-components) "Displ_1" and "Displ_2" in form of three-dimensional
% arrays. In input is also required the measurement grid "X","Y" and the
% acquisition frequency "facq" and SSF.
% The dimensions of "Displ_1" are (J,I,Z), where J and I are the number of
% rows and columns.  The third dimension is used to distinguish between
% the x- and y-component of the displacement.
% The routine output is a four-dimensional array containing the time
% super-sampled data serries "Displ_TSS" of dimensions (J,I,SSF,2).

%function time_super_sample(Displx,Disply,facq,SSF,grid_scale)

    count = 0;
    clear Displx_new Disply_new
    [X,Y] = meshgrid(1:size(Disply,2),1:size(Displx,1));
    [Xn,Yn] = meshgrid(1:1/grid_scale:size(Disply,2),1:1/grid_scale:size(Displx,1));
    X = grid_scale*X; Y = grid_scale*Y;
    
    % time separation between measurements
    DT = 1/facq;
    
    % time separation after super-sampling
    dt = 1/facq/SSF;
%%
% for i =1:size(Displx,3)
for i = 1:size(Displx,3)-1
   
    Displx_1 = Displx(:,:,i);
    Displx_2 = Displx(:,:,i+1);
    
    Disply_1 = Disply(:,:,i);
    Disply_2 = Disply(:,:,i+1);
    
    % initialize the super-sampled sequence
    Displx_TSS = zeros(size(X,1),size(X,2),SSF);
    Disply_TSS = zeros(size(Y,1),size(Y,2),SSF);
    
    for n = 0:SSF-1
        % first approximation

        % determine the positions of the fluid parcel at a subsequent time
        % instant based on the convective velocity field estimated from first
        % field
        Xdt1 = X-Displx_1*(DT/dt)*(n/SSF);
        Ydt1 = Y-Disply_1*(DT/dt)*(n/SSF);
        
        % determine the positions of the fluid parcel at a subsequent time
        % instant based on the convective velocity field estimated from first
        % field
        Xdt2 = X+Displx_2*(DT/dt)*(1-n/SSF);
        Ydt2 = Y+Disply_2*(DT/dt)*(1-n/SSF);
        
        % evaluate the velocity at the arrival position from spatial
        % interpolation
        disp('Calculating first approx. on U')
        Udt1 = griddata(X,Y,Displx_1,Xdt1,Ydt1,'v4');
        Udt2 = griddata(X,Y,Displx_2,Xdt2,Ydt2,'v4');
        
        disp('Calculating first approx. on V')
        Vdt1 = griddata(X,Y,Disply_1,Xdt1,Ydt1,'v4');
        Vdt2 = griddata(X,Y,Disply_2,Xdt2,Ydt2,'v4');
        
        % weighted average of velocity from backward and forward evaluation
        Utss = Udt1*(1-n/SSF)+Udt2*(n/SSF);
        Vtss = Vdt1*(1-n/SSF)+Vdt2*(n/SSF);
        
        % second approximation
        % improved estimate of convective velocity
        Xdt1 = X-(Displx_1+Udt1)/2*(DT/dt)*(n/SSF);
        Ydt1 = Y-(Disply_1+Vdt1)/2*(DT/dt)*(n/SSF);
        
        Xdt2 = X+(Displx_2+Udt2)/2*(DT/dt)*(1-n/SSF);
        Ydt2 = Y+(Disply_2+Vdt2)/2*(DT/dt)*(1-n/SSF);
        
        disp('Calculating second approx. on U')
        Udt1 = griddata(X,Y,Displx_1,Xdt1,Ydt1,'v4');
        Udt2 = griddata(X,Y,Displx_2,Xdt2,Ydt2,'v4');
        
        disp('Calculating second approx. on V')
        Vdt1 = griddata(X,Y,Disply_1,Xdt1,Ydt1,'v4');
        Vdt2 = griddata(X,Y,Disply_2,Xdt2,Ydt2,'v4');
        
        Utss = Udt1*(1-n/SSF)+Udt2*(n/SSF);
        Vtss = Vdt1*(1-n/SSF)+Vdt2*(n/SSF);
        
        % store result into sequence array
        Displx_TSS(:,:,n+1)=(Utss);
        Disply_TSS(:,:,n+1)=(Vtss);
    end;
    
    for t = 1:size(Displx_TSS,3)
        count = count+1
        Displx_new(:,:,count) = Displx_TSS(:,:,t);
        Disply_new(:,:,count) = Disply_TSS(:,:,t);
    end;
end;