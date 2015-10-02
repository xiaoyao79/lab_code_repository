function [Output] = ThreePIcode(Data)
%function [Output]=ThreePIcode(Data)
% 3Pi function to make the flow calculations.  flowvarsteady [Jmax x Imax x numbervariables] is the time average quantities of
% flowvarunsteady [Jmax x Imax x numframes x numbervariables], some
% variables such as skewness and kurtosis are only defined for steady
% calculations.  to add new variables, go the the end of this file and append a new element along the
% fourth dimension to flowvarunsteady and append along the third dim for flowvarsteady.
% Also modify (add new element to vector) stdyvarmat and unstdyvarmat with the
% number of new quantities (this is not necessary, it is for memory
% preallocation).  also, when adding new variables, add the appropriate
% string labels to varcellstdy and varcellunstdy, they are printed in the
% header output files
% note, the calcultaions take place assuming x varies with
% columns and y varies with the rows.
% varout are the extra variables from the plt file that were loaded
% (variables other than x,y,u,v)
% um,vm,wm are the mean velocity values if you are calling this function
% one frame at a time, they should be empty if you are calculating multiple
% frames.

% Last Updated - 05/07/09 Sam Raben
%              -
To = clock;
% clc
if nargin == 0
    return
end
Output = Data;

% if Data.calc.gradientsvar == 0 && Data.calc.vortexidvar == 1
%     fprintf('Vortex ID can not be preformed with out Gradient Methods Selected\n')
% end

% define the flow properties
charvel    = str2double(Data.phys.charvel);          % Characteristic Velocity
charlen    = str2double(Data.phys.charlen);          % Characteristic Length
freq       = str2double(Data.phys.freq);             % Sample Frequence
micrperpix = str2double(Data.phys.micrperpix);       % Resolution (um/pix)
pulsesep   = str2double(Data.phys.pulsesep);         % Pulse Separation (us)
kinvisc    = str2double(Data.phys.kinvisc);          % Kinimatic Viscosity (m^2/2)
xoffset    = Data.phys.xoffset;          %
yoffset    = Data.phys.yoffset;          %
zoffset    = Data.phys.zoffset;          %#ok
unittype   = Data.phys.unitvar;          % 1 = pixel units, 2 = m/s

if ~Data.calc.vorticityvar && ~Data.calc.gradientsvar
    Data.calc.gradientsvar = 1;
    fprintf('Turning on Gradients to calculate Vorticity\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
casename   = Data.plt_info.pltname;
dirloc     = Data.plt_info.pltdir;
startframe = str2double(Data.plt_info.startframe);
endframe   = str2double(Data.plt_info.endframe);
numzeros   = str2double(Data.plt_info.numzeros);
framestep  = str2double(Data.plt_info.skipframes);
%varlist    = Data.calculations.varlist;
if isfield(Data.plt_info,'loadpod')
    if Data.plt_info.loadpod == 1

        [PLT.x,PLT.y,PLT.u,PLT.v,PLT.extravars,PLT.varlistnew]=read_pltmod(casename,dirloc,startframe,startframe,numzeros,framestep);
        POD = load(Data.plt_info.podfile);
        PLT.umn = POD.u_mean;
        PLT.vmn = POD.v_mean;
        U_SIZE = size(POD.V1);
        recnum = find((cumsum(POD.D)/sum(POD.D))>Data.plt_info.podReconPer,1,'first');
        PLT.u   = reshape(reshape(POD.V1(:,:,1:recnum),U_SIZE(1)*U_SIZE(2),recnum)*POD.A(1:recnum,1:U_SIZE(3)),U_SIZE(1),U_SIZE(2),U_SIZE(3))+repmat(POD.u_mean,[1 1 U_SIZE(3)]);
        PLT.v   = reshape(reshape(POD.V2(:,:,1:recnum),U_SIZE(1)*U_SIZE(2),recnum)*POD.A(1:recnum,1:U_SIZE(3)),U_SIZE(1),U_SIZE(2),U_SIZE(3))+repmat(POD.v_mean,[1 1 U_SIZE(3)]);
        clear POD
        PLT.varlist = {};
    else
        %[PLT.x,PLT.y,PLT.u,PLT.v,PLT.extravars,PLT.varlistnew]=read_pltmod(casename,dirloc,startframe,endframe,numzeros,framestep,varlist);%#ok
        [PLT.x,PLT.y,PLT.u,PLT.v,PLT.extravars,PLT.varlistnew]=read_pltmod(casename,dirloc,startframe,endframe,numzeros,framestep);
        PLT.umn = mean(PLT.u,3);
        PLT.vmn = mean(PLT.v,3);
    end
else
    %[PLT.x,PLT.y,PLT.u,PLT.v,PLT.extravars,PLT.varlistnew]=read_pltmod(casename,dirloc,startframe,endframe,numzeros,framestep,varlist);%#ok
    [PLT.x,PLT.y,PLT.u,PLT.v,PLT.extravars,PLT.varlistnew]=read_pltmod(casename,dirloc,startframe,endframe,numzeros,framestep);
    PLT.umn = mean(PLT.u,3);
    PLT.vmn = mean(PLT.v,3);
end

T1 = etime(clock,To);
fprintf('Load Time...                    %0.2i:%0.2i.%0.0f\n',floor(T1/60),floor(rem(T1,60)),...
    rem(T1,60)-floor(rem(T1,60)))
To = clock;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparing Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Jdat,Idat,numframes]=size(PLT.u);

wflag=0;
for j=1:length(PLT.varlistnew)
    if strcmpi(PLT.varlistnew{j},'w')         % determine if the third velocity component data is present
        PLT.w     = PLT.extravars{j};
        PLT.u     = PLT.v;
        PLT.v     = PLT.extravars{j-1};
        PLT.wmn   = mean(PLT.w,3);
        PLT.umn   = mean(PLT.u,3);
        PLT.vmn   = mean(PLT.v,3);
        wflag=1;
    else
        
    end
    
end
if ~wflag
    PLT.w = zeros(size(PLT.u,1),size(PLT.u,2),1);
    PLT.wmn = zeros(size(PLT.w));
end
    

time=(0:numframes-1)/freq;      % define the time at each frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert in Physical Units from pixels if Necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if unittype==1          % normalize for input units of pixels and delta t
    PLT.x       = (PLT.x*micrperpix*1e-6+xoffset);
    PLT.y       = (PLT.y*micrperpix*1e-6+yoffset);
    PLT.u       = PLT.u/pulsesep*micrperpix;
    PLT.v       = PLT.v/pulsesep*micrperpix;
%    kinvisc     = kinvisc;
%    time        = time;
    PLT.umn     = (PLT.umn)*(1/pulsesep)*micrperpix;
    PLT.vmn     = (PLT.vmn)*(1/pulsesep)*micrperpix;
    if wflag>0
        PLT.w       = PLT.w/pulsesep*micrperpix;
        PLT.wmn     = (PLT.wmn/pulsesep*micrperpix);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% determine delta x and delta y spacing, and dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=PLT.x(1,2)-PLT.x(1,1);
dy=PLT.y(2,1)-PLT.y(1,1);

uprime = zeros(size(PLT.u));
vprime = zeros(size(PLT.v));
if wflag > 0
    wprime = zeros(size(PLT.w));
end

% initialize for speed
for ipm = 1:numel(PLT.u(1,1,:))
    uprime(:,:,ipm) = PLT.u(:,:,ipm) - PLT.umn;
    vprime(:,:,ipm) = PLT.v(:,:,ipm) - PLT.vmn;
    if wflag > 0
        wprime(:,:,ipm)=PLT.w(:,:,ipm) - PLT.wmn;
    end
end

numvarssteady = Data.calc.velmagvar + Data.calc.vorticityvar + Data.calc.strainratevar + Data.calc.turbKEvar + ...
    Data.calc.totalKEvar + Data.calc.reynoldsstressvar + (Data.calc.turbstats*5) + Data.calc.diss1term + ...
    Data.calc.diss5term + Data.calc.dissLES + Data.calc.dissCOMP + (Data.calc.gradientsvar * 8);

numvarsunsteady = Data.calc.velmagvar + Data.calc.vorticityvar + Data.calc.strainratevar + Data.calc.turbKEvar + ...
    Data.calc.totalKEvar  + Data.calc.diss1term + ...
    Data.calc.diss5term + Data.calc.dissLES + Data.calc.dissCOMP + Data.calc.gradientsvar * 5;

if Data.calc.steady
    Output.flowvarsteady   = zeros(Jdat,Idat,1,numvarssteady,'single');
end

if Data.calc.unsteady == 1 && (Data.output.fileformat == 1 || Data.output.fileformat == 3)
    Output.flowvarunsteady = zeros(Jdat,Idat,numframes,numvarsunsteady,'single');
end

singlematdat = (Data.output.fileformat == 1 || Data.output.fileformat == 3 || Data.output.fileformat == 4);

if Data.calc.gradientsvar
    fprintf('Calculating Gradients...')
    if Data.calc.unsteady && Data.calc.gradientsvar
        Data.calc.gradienttype = 1;
    end

    GMETHOD = {'order2','rich4','chapra'};

    if Data.calc.steady && ~Data.calc.unsteady %~Data.calc.gradienttype
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gradient of the Average
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Grad.dudx=zeros(Jdat,Idat,'single');
        Grad.dvdx=zeros(Jdat,Idat,'single');
        Grad.dudy=zeros(Jdat,Idat,'single');
        Grad.dvdy=zeros(Jdat,Idat,'single');
        Grad.dudz=zeros(Jdat,Idat,'single');
        Grad.dvdz=zeros(Jdat,Idat,'single');
        Grad.dwdx=zeros(Jdat,Idat,'single');
        Grad.dwdy=zeros(Jdat,Idat,'single');


        switch lower(GMETHOD{Data.calc.gradientmethod})
            case'order2'
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_order2(PLT.umn,PLT.vmn,PLT.wmn,dx,dy);
            case 'rich4'
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_rich4(PLT.umn,PLT.vmn,PLT.wmn,dx,dy);
            case 'chapra'
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_chapra(PLT.umn,PLT.vmn,PLT.wmn,dx,dy);
        end

    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Average of the Gradient
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Grad.dudx=zeros(Jdat,Idat,numframes,'single');
        Grad.dvdx=zeros(Jdat,Idat,numframes,'single');
        Grad.dudy=zeros(Jdat,Idat,numframes,'single');
        Grad.dvdy=zeros(Jdat,Idat,numframes,'single');
        Grad.dudz=zeros(Jdat,Idat,numframes,'single');
        Grad.dvdz=zeros(Jdat,Idat,numframes,'single');
        Grad.dwdx=zeros(Jdat,Idat,numframes,'single');
        Grad.dwdy=zeros(Jdat,Idat,numframes,'single');

        for j=1:numframes
            um=PLT.u(:,:,j);
            vm=PLT.v(:,:,j);
            if wflag > 0
                wm=PLT.w(:,:,j);
            else
                wm=zeros(size(PLT.u(:,:,1)));
            end

            % calculate gradients using chosen method
            switch lower(GMETHOD{Data.calc.gradientmethod})
                case'order2'
                    [Grad.dudx(:,:,j),Grad.dudy(:,:,j),Grad.dvdx(:,:,j),Grad.dvdy(:,:,j),Grad.dwdx(:,:,j),Grad.dwdy(:,:,j)] =...
                        gradient_dudxdvdy_order2(um,vm,wm,dx,dy);
                case 'rich4'
                    [Grad.dudx(:,:,j),Grad.dudy(:,:,j),Grad.dvdx(:,:,j),Grad.dvdy(:,:,j),Grad.dwdx(:,:,j),Grad.dwdy(:,:,j)] =...
                        gradient_dudxdvdy_rich4(um,vm,wm,dx,dy);
                case 'chapra'
                    [Grad.dudx(:,:,j),Grad.dudy(:,:,j),Grad.dvdx(:,:,j),Grad.dvdy(:,:,j),Grad.dwdx(:,:,j),Grad.dwdy(:,:,j)] =...
                        gradient_dudxdvdy_chapra(um,vm,wm,dx,dy);
            end
        end
    end
end


%%
%%%%%% put x,y,u,v and extra variables into the 4D arrays
if Data.calc.steady
    if Data.phys.nondim
        Output.flowvarsteady(:,:,:,1)   = PLT.x./charlen;
        Output.flowvarsteady(:,:,:,2)   = PLT.y./charlen;
        Output.flowvarsteady(:,:,:,3)   = PLT.umn./charvel;
        Output.flowvarsteady(:,:,:,4)   = PLT.vmn./charvel;
        Output.varcellstdy={'X/L' 'Y/L' 'U/Uo' 'V/Uo'};
    else
        Output.flowvarsteady(:,:,:,1)   = PLT.x;
        Output.flowvarsteady(:,:,:,2)   = PLT.y;
        Output.flowvarsteady(:,:,:,3)   = PLT.umn;
        Output.flowvarsteady(:,:,:,4)   = PLT.vmn;
        Output.varcellstdy={'X' 'Y' 'U' 'V'};
    end
end

if Data.calc.unsteady && singlematdat
    if Data.phys.nondim
        Output.flowvarunsteady(:,:,:,1) = repmat(PLT.x,[1 1 numframes])./charlen;
        Output.flowvarunsteady(:,:,:,2) = repmat(PLT.y,[1 1 numframes])./charlen;
        Output.flowvarunsteady(:,:,:,3) = PLT.u./charvel;
        Output.flowvarunsteady(:,:,:,4) = PLT.v./charvel;
        Output.varcellunstdy={'X/L' 'Y/L' 'U/Uo' 'V/Uo'};
    else
        Output.flowvarunsteady(:,:,:,1) = repmat(PLT.x,[1 1 numframes]);
        Output.flowvarunsteady(:,:,:,2) = repmat(PLT.y,[1 1 numframes]);
        Output.flowvarunsteady(:,:,:,3) = PLT.u;
        Output.flowvarunsteady(:,:,:,4) = PLT.v;
        Output.varcellunstdy={'X' 'Y' 'U' 'V'};
    end
end

k=5;        % inititalize counters
if Data.calc.steady
    for j=1:length(PLT.varlistnew)     % put in the extra loaded variables from plt
        if strcmpi(PLT.varlistnew{j},'w')%wflag > 0                 % determine if one of the extra variables is w
            if Data.phys.nondim
                Output.flowvarsteady(:,:,:,k)   = PLT.wmn./charvel;
                Output.varcellstdy{k}           = 'W/Uo';
            else
                Output.flowvarsteady(:,:,:,k)   = PLT.wmn;
                Output.varcellstdy{k}           = 'W';
            end
            k=k+1;
        elseif ~strcmpi(PLT.varlistnew{j},'u') && ~strcmpi(PLT.varlistnew{j},'v')
            Output.flowvarsteady(:,:,:,k)   = mean(PLT.extravars{1},3);
            Output.varcellstdy{k}           = PLT.varlistnew{j};
            k=k+1;
        end
    end
end

c=5;       % inititalize counters
if Data.calc.unsteady && singlematdat
    for j=1:length(PLT.varlistnew)     % put in the extra loaded variables from plt
        if strcmpi(PLT.varlistnew{j},'w')                % determine if one of the extra variables is w
            if Data.phys.nondim
                Output.flowvarunsteady(:,:,:,c) = PLT.w./charvel;
                Output.varcellunstdy{c}         = 'W/Uo';
            else
                Output.flowvarunsteady(:,:,:,c) = PLT.w;
                Output.varcellunstdy{c}         = 'W';
            end
            c=c+1;
        elseif ~strcmpi(PLT.varlistnew{j},'u') && ~strcmpi(PLT.varlistnew{j},'v') && ~strcmpi(PLT.varlistnew{j},'v') ...
            && ~strcmpi(PLT.varlistnew{j},'U1') && ~strcmpi(PLT.varlistnew{j},'V1') && ~strcmpi(PLT.varlistnew{j},'U2') ...
            && ~strcmpi(PLT.varlistnew{j},'V2') && ~strcmpi(PLT.varlistnew{j},'U3') && ~strcmpi(PLT.varlistnew{j},'V3') ...
            && ~strcmpi(PLT.varlistnew{j},'C1') && ~strcmpi(PLT.varlistnew{j},'C2') && ~strcmpi(PLT.varlistnew{j},'C3') ...
            && ~strcmpi(PLT.varlistnew{j},'D1') && ~strcmpi(PLT.varlistnew{j},'D2') && ~strcmpi(PLT.varlistnew{j},'D3')
            try
%             Output.flowvarunsteady(:,:,:,c) = PLT.extravars{j};
%             Output.varcellunstdy{c}         = PLT.varlistnew{j};
%             c=c+1;
            catch ME
                fprintf('%s\n',ME.message)
                keyboard
            end
        end
    end
end

k = k-1;c=c-1;
fprintf('Done')
T1 = etime(clock,To);
fprintf('    %0.2i:%0.2i.%0.0f\n',floor(T1/60),floor(rem(T1,60)),...
    rem(T1,60)-floor(rem(T1,60)))
To = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calculating Statistics...')
%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity Magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.velmagvar
    if wflag == 0
        velmag=sqrt(PLT.u.^2+PLT.v.^2);
    else
        velmag=sqrt(PLT.u.^2+PLT.v.^2+PLT.w.^2);
    end

    c=c+1;k=k+1;
    if Data.calc.steady
        if Data.phys.nondim
            Output.flowvarsteady(:,:,:,k)=mean(velmag,3)./charvel;
            Output.varcellstdy{k}='Velmag/Uo';
        else
            Output.flowvarsteady(:,:,:,k)=mean(velmag,3);
            Output.varcellstdy{k}='Velmag';
        end        
    end

    if Data.calc.unsteady && singlematdat
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=velmag./charvel;
            Output.varcellunstdy{c}='Velmag/Uo';
        else
            Output.flowvarunsteady(:,:,:,c)=velmag;
            Output.varcellunstdy{c}='Velmag';
        end
    end
    clear velmag;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Vorticity
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.vorticityvar
    vort=Grad.dvdx-Grad.dudy;
    c=c+1;k=k+1;
    if Data.calc.steady
        if Data.phys.nondim
            Output.flowvarsteady(:,:,:,k)=mean(vort,3).*charlen/charvel;
            Output.varcellstdy{k}='Vorticity*L/Uo';
        else
            Output.flowvarsteady(:,:,:,k)=mean(vort,3);
            Output.varcellstdy{k}='Vorticity';
        end
    end
    if Data.calc.unsteady && singlematdat
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=vort.*charlen/charvel;
            Output.varcellunstdy{c}='Vorticity*L/Uo';
        else
            Output.flowvarunsteady(:,:,:,c)=vort;
            Output.varcellunstdy{c}='Vorticity';
        end
    end
    clear vort;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Strain rate, off diagonal term of 2x2 tensor
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.strainratevar
    strainrate=Grad.dudy+Grad.dvdx;
    c=c+1;k=k+1;
    if Data.calc.steady
        if Data.phys.nondim
            Output.flowvarsteady(:,:,:,k)=mean(strainrate,3).*charlen/charvel;
            Output.varcellstdy{k}='Strain Rate*L/Uo';
        else
            Output.flowvarsteady(:,:,:,k)=mean(strainrate,3);
            Output.varcellstdy{k}='Strain Rate';
        end
    end
    if Data.calc.unsteady && singlematdat
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=strainrate.*charlen/charvel;
            Output.varcellunstdy{c}='Strain Rate*L/Uo';
        else
            Output.flowvarunsteady(:,:,:,c)=strainrate;
            Output.varcellunstdy{c}='Strain Rate';
        end
    end
    clear strainrate;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Gradients
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.gradientsvar        % gradients
    c=c+4;k=k+4;
    if Data.calc.steady
        if Data.phys.nondim
            Output.flowvarsteady(:,:,:,k-3)=mean(Grad.dudx,3).*charlen/charvel;
            Output.flowvarsteady(:,:,:,k-2)=mean(Grad.dvdx,3).*charlen/charvel;
            Output.flowvarsteady(:,:,:,k-1)=mean(Grad.dudy,3).*charlen/charvel;
            Output.flowvarsteady(:,:,:,k)=mean(Grad.dvdy,3).*charlen/charvel;
            Output.varcellstdy{k-3}='dudx*L/Uo';
            Output.varcellstdy{k-2}='dvdx*L/Uo';
            Output.varcellstdy{k-1}='dudy*L/Uo';
            Output.varcellstdy{k}='dvdy*L/Uo';
            if wflag>0
                Output.flowvarsteady(:,:,:,k+1)=mean(Grad.dwdx,3).*charlen/charvel;
                Output.flowvarsteady(:,:,:,k+2)=mean(Grad.dwdy,3).*charlen/charvel;
                Output.varcellstdy{k+1}='dwdx*L/Uo';
                Output.varcellstdy{k+2}='dwdy*L/Uo';
                k = k+2;
            end
        else
            Output.flowvarsteady(:,:,:,k-3)=mean(Grad.dudx,3);
            Output.flowvarsteady(:,:,:,k-2)=mean(Grad.dvdx,3);
            Output.flowvarsteady(:,:,:,k-1)=mean(Grad.dudy,3);
            Output.flowvarsteady(:,:,:,k)=mean(Grad.dvdy,3);
            Output.varcellstdy{k-3}='dudx';
            Output.varcellstdy{k-2}='dvdx';
            Output.varcellstdy{k-1}='dudy';
            Output.varcellstdy{k}='dvdy';
            if wflag>0
                Output.flowvarsteady(:,:,:,k+1)=mean(Grad.dwdx,3);
                Output.flowvarsteady(:,:,:,k+2)=mean(Grad.dwdy,3);
                Output.varcellstdy{k+1}='dwdx';
                Output.varcellstdy{k+2}='dwdy';
                k = k+2;
            end
        end
    end

    if Data.calc.unsteady && singlematdat
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c-3)=Grad.dudx.*charlen/charvel;
            Output.flowvarunsteady(:,:,:,c-2)=Grad.dvdx.*charlen/charvel;
            Output.flowvarunsteady(:,:,:,c-1)=Grad.dudy.*charlen/charvel;
            Output.flowvarunsteady(:,:,:,c)=Grad.dvdy.*charlen/charvel;
            Output.varcellunstdy{c-3}='dudx*L/Uo';
            Output.varcellunstdy{c-2}='dvdx*L/Uo';
            Output.varcellunstdy{c-1}='dudy*L/Uo';
            Output.varcellunstdy{c}='dvdy*L/Uo';
            if wflag > 0;
                Output.flowvarunsteady(:,:,:,c+1)=Grad.dwdx.*charlen/charvel;
                Output.flowvarunsteady(:,:,:,c+2)=Grad.dwdy.*charlen/charvel;
                Output.varcellunstdy{c+1}='dwdx*L/Uo';
                Output.varcellunstdy{c+2}='dwdy*L/Uo';
                c = c+2;
            end
        else
            Output.flowvarunsteady(:,:,:,c-3)=Grad.dudx;
            Output.flowvarunsteady(:,:,:,c-2)=Grad.dvdx;
            Output.flowvarunsteady(:,:,:,c-1)=Grad.dudy;
            Output.flowvarunsteady(:,:,:,c)=Grad.dvdy;
            Output.varcellunstdy{c-3}='dudx';
            Output.varcellunstdy{c-2}='dvdx';
            Output.varcellunstdy{c-1}='dudy';
            Output.varcellunstdy{c}='dvdy';
            if wflag > 0;
                Output.flowvarunsteady(:,:,:,c+1)=Grad.dwdx;
                Output.flowvarunsteady(:,:,:,c+2)=Grad.dwdy;
                Output.varcellunstdy{c+1}='dwdx';
                Output.varcellunstdy{c+2}='dwdy';
                c = c+2;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Turbulent Kinetic Energy
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.turbKEvar 
    if wflag == 0;
        turbke = uprime.^2+vprime.^2;
    else
        turbke = uprime.^2+vprime.^2+wprime.^2;
    end

    c=c+1;k=k+1;
    if Data.calc.steady
        if Data.phys.nondim
            Output.flowvarsteady(:,:,:,k)=mean(turbke,3)./(charvel^2);
            Output.varcellstdy{k}='TurbulentKE/Uo^2';
        else
            Output.flowvarsteady(:,:,:,k)=mean(turbke,3);
            Output.varcellstdy{k}='TurbulentKE';
        end
    end
    
    if Data.calc.unsteady && singlematdat
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=turbke./(charvel^2);
            Output.varcellunstdy{c}='TurbulentKE/Uo^2';
        else
            Output.flowvarunsteady(:,:,:,c)=turbke;
            Output.varcellunstdy{c}='TurbulentKE';
        end
    end
    clear turbke;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Kinetic Energy
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.totalKEvar 
    if wflag == 0;
        totalke=PLT.u.^2+PLT.v.^2;
    else
        totalke=PLT.u.^2+PLT.v.^2+PLT.w.^2;
    end

    c=c+1;k=k+1;
    if Data.calc.steady
        if Data.phys.nondim
            Output.flowvarsteady(:,:,:,k)=mean(totalke,3)./(charvel^2);
            Output.varcellstdy{k}='TotalKE/Uo^2';
        else
            Output.flowvarsteady(:,:,:,k)=mean(totalke,3);
            Output.varcellstdy{k}='TotalKE';
        end
    end

    if Data.calc.unsteady && singlematdat
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=totalke./(charvel^2);
            Output.varcellunstdy{c}='TotalKE/Uo^2';
        else
            Output.flowvarunsteady(:,:,:,c)=totalke;
            Output.varcellunstdy{c}='TotalKE';
        end
    end
    clear totalke;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Reynolds Stresses
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.reynoldsstressvar
        c=c+3;k=k+3;
        if Data.calc.steady
            if Data.phys.nondim
                Output.flowvarsteady(:,:,:,k-2) = mean(uprime.^2,3)./(charvel^2);
                Output.flowvarsteady(:,:,:,k-1) = mean(vprime.^2,3)./(charvel^2);
                Output.flowvarsteady(:,:,:,k)   = mean(uprime.*vprime,3)./(charvel^2);
                Output.varcellstdy{k-2}         = 'u''^2/Uo^2';
                Output.varcellstdy{k-1}         = 'v''^2/Uo^2';
                Output.varcellstdy{k}           = 'u''v''/Uo^2';
                if wflag > 0
                    Output.flowvarsteady(:,:,:,k+1)=mean(wprime.^2,3)./(charvel^2);
                    Output.varcellstdy{k+1}='w''^2/Uo^2';
                    k=k+1;
                end
            else
                Output.flowvarsteady(:,:,:,k-2) = mean(uprime.^2,3);
                Output.flowvarsteady(:,:,:,k-1) = mean(vprime.^2,3);
                Output.flowvarsteady(:,:,:,k)   = mean(uprime.*vprime,3);
                Output.varcellstdy{k-2}         = 'u''^2';
                Output.varcellstdy{k-1}         = 'v''^2';
                Output.varcellstdy{k}           = 'u''v''';
                if wflag > 0
                    Output.flowvarsteady(:,:,:,k+1)=mean(wprime.^2,3);
                    Output.varcellstdy{k+1}='w''^2';
                    k=k+1;
                end
            end
        end

        if Data.calc.unsteady && singlematdat
            if Data.phys.nondim
                Output.flowvarunsteady(:,:,:,c-2)=(uprime.^2)./(charvel^2);
                Output.flowvarunsteady(:,:,:,c-1)=(vprime.^2)./(charvel^2);
                Output.flowvarunsteady(:,:,:,c)=(uprime.*vprime)./(charvel^2);
                Output.varcellunstdy{c-2}='u''^2/Uo^2';
                Output.varcellunstdy{c-1}='v''^2/Uo^2';
                Output.varcellunstdy{c}='u''v''/Uo^2';
                if wflag > 0
                    Output.flowvarunsteady(:,:,:,c+1)=(wprime.^2)./(charvel^2);
                    Output.varcellunstdy{c+1}='w''^2/Uo^2';
                    c=c+1;
                end
            else
                Output.flowvarunsteady(:,:,:,c-2)=uprime.^2;
                Output.flowvarunsteady(:,:,:,c-1)=vprime.^2;
                Output.flowvarunsteady(:,:,:,c)=uprime.*vprime;
                Output.varcellunstdy{c-2}='u''^2';
                Output.varcellunstdy{c-1}='v''^2';
                Output.varcellunstdy{c}='u''v''';
                if wflag > 0
                    Output.flowvarunsteady(:,:,:,c+1)=wprime.^2;
                    Output.varcellunstdy{c+1}='w''^2';
                    c=c+1;
                end                
            end
        end
    clear reystressuprime reystressvprime reystresswprime reystressupvp reystressupwp reystressvpwp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Turbulent Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.turbstats
    if Data.calc.steady
        if ~isempty(PLT.umn)
            urms = std(PLT.u,1,3);
            vrms = std(PLT.v,1,3);
            uskew = skewness(PLT.u,1,3);
            vskew = skewness(PLT.v,1,3);
            ukurt = kurtosis(PLT.u,1,3);
            vkurt = kurtosis(PLT.v,1,3);
            if wflag >0
                wrms = std(PLT.w,1,3);
                wskew = skewness(PLT.w,1,3);
                wkurt = kurtosis(PLT.w,1,3);
            end
        else
            urms=uprime.^2;
            vrms=vprime.^2;
            uskew=uprime.^3;
            vskew=vprime.^3;
            ukurt=mean(uprime.^4,3);
            vkurt=mean(vprime.^4,3);
        end

        % remove the NaN's, these arise from taking the skewness
        % or kurtosis of a vector of all zeros.
        [rrr,ccc]=find(isnan(uskew)==1 | isnan(vskew)==1);
        uskew(rrr,ccc)=0;
        vskew(rrr,ccc)=0;
        ukurt(rrr,ccc)=0;
        vkurt(rrr,ccc)=0;
        if wflag
            wskew(isnan(wskew)) = 0;
            wkurt(isnan(wkurt)) = 0;
        end
        clear rrr ccc

        if wflag == 0;
            k=k+6;
            if Data.phys.nondim
                Output.flowvarsteady(:,:,:,k-5)=urms./charvel;
                Output.flowvarsteady(:,:,:,k-4)=vrms./charvel;
                Output.varcellstdy{k-5}='Urms/Uo';
                Output.varcellstdy{k-4}='Vrms/Uo';
            else
                Output.flowvarsteady(:,:,:,k-5)=urms;
                Output.flowvarsteady(:,:,:,k-4)=vrms;
                Output.varcellstdy{k-5}='Urms';
                Output.varcellstdy{k-4}='Vrms';
            end
            
            Output.flowvarsteady(:,:,:,k-3)=uskew;
            Output.flowvarsteady(:,:,:,k-2)=vskew;
            Output.flowvarsteady(:,:,:,k-1)=ukurt;
            Output.flowvarsteady(:,:,:,k)=vkurt;
            Output.varcellstdy{k-3}='Uskew';
            Output.varcellstdy{k-2}='Vskew';
            Output.varcellstdy{k-1}='Ukurt';
            Output.varcellstdy{k}='Vkurt';
        else
            k = k+9;
            if Data.phys.nondim
                Output.flowvarsteady(:,:,:,k-8)= urms./charvel;
                Output.flowvarsteady(:,:,:,k-7)= vrms./charvel;
                Output.flowvarsteady(:,:,:,k-6)= wrms./charvel;
                Output.varcellstdy{k-8}='Urms/Uo';
                Output.varcellstdy{k-7}='Vrms/Uo';
                Output.varcellstdy{k-6}='Wrms/Uo';
            else
                Output.flowvarsteady(:,:,:,k-8)= urms;
                Output.flowvarsteady(:,:,:,k-7)= vrms;
                Output.flowvarsteady(:,:,:,k-6)= wrms;
                Output.varcellstdy{k-8}='Urms';
                Output.varcellstdy{k-7}='Vrms';
                Output.varcellstdy{k-6}='Wrms';                
            end
            
            Output.flowvarsteady(:,:,:,k-5)= uskew;
            Output.flowvarsteady(:,:,:,k-4)= vskew;
            Output.flowvarsteady(:,:,:,k-3)= wskew;
            Output.flowvarsteady(:,:,:,k-2)= ukurt;
            Output.flowvarsteady(:,:,:,k-1)= vkurt;
            Output.flowvarsteady(:,:,:,k)  = wkurt;
            Output.varcellstdy{k-5}='Uskew';
            Output.varcellstdy{k-4}='Vskew';
            Output.varcellstdy{k-3}='Wskew';
            Output.varcellstdy{k-2}='Ukurt';
            Output.varcellstdy{k-1}='Vkurt';
            Output.varcellstdy{k}='Wkurt';
        end
        %     end
        clear urms vrms wrms uskew vskew wskew ukurt vkurt wkurt;
    end
end
fprintf('Done')
T1 = etime(clock,To);
fprintf('   %0.2i:%0.2i.%0.0f\n',floor(T1/60),floor(rem(T1,60)),...
    rem(T1,60)-floor(rem(T1,60)))
To = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dissipation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 Term Dissipation
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.diss1term+Data.calc.diss5term+Data.calc.dissLES + Data.calc.dissCOMP > 0
    fprintf('Calculating Dissipations...')
end
if Data.calc.diss1term
    if Data.calc.diss1var == 1
        term1=mean(Grad.dudx.^2,3);
    elseif Data.calc.diss1var == 2;
        term1=mean(Grad.dvdy.^2,3);
    else
        term1=mean(Grad.dwdz.^2,3);
    end

    k=k+1;
    if Data.phys.nondim
        Output.flowvarsteady(:,:,:,k) = (kinvisc*15*term1)./(charlen/(charvel ^3));
        Output.varcellstdy{k}='dissipation(1T)*L/Uo^3';
    else
        Output.flowvarsteady(:,:,:,k) = kinvisc*15*term1;
        Output.varcellstdy{k}='dissipation(1T)';
    end
    clear term1 dissipgrad
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 Term Dissipation
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.diss5term
    term1 = mean(Grad.dudx.^2,3);
    term2 = mean(Grad.dvdy.^2,3);
    term3 = mean(Grad.dudy.^2,3);
    term4 = mean(Grad.dvdx.^2,3);
    term5 = mean(Grad.dudy.*Grad.dvdx,3);

    dissipgrad = kinvisc*(2.*term1+2.*term2+3.*term3+3.*term4+2.*term5);
    
    k=k+1;
    if Data.phys.nondim
        Output.flowvarsteady(:,:,:,k)=dissipgrad./(charlen/(charvel^3));
        Output.varcellstdy{k}='dissipation(5T)*L/Uo^3';
    else
        Output.flowvarsteady(:,:,:,k)=dissipgrad;
        Output.varcellstdy{k}='dissipation(5T)';
    end
    clear term1 term2 term3 term4 term5 dissipgrad;
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% LES Term Dissipation
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.dissLES
    %fprintf('Calculating LES Dissipation...')
    coef = str2double(Data.calc.dissLEScoef);
    filt = str2double(Data.calc.dissLESfilt);

    %preallocate main TR and TA dissipation arrays
    epsilon_TR=zeros(size(Grad.dudx));

    %main loop for calculating dissipation at each time step
    for i=1:size(Grad.dudx,3)
        %directly calculate the strain rate tensor components from the time-
        %resolved two-dimensional gradient fields
        S11=0.5*(Grad.dudx(:,:,i)+Grad.dudx(:,:,i));
        S12=0.5*(Grad.dudy(:,:,i)+Grad.dvdx(:,:,i));
        S21=0.5*(Grad.dvdx(:,:,i)+Grad.dudy(:,:,i));
        S22=0.5*(Grad.dvdy(:,:,i)+Grad.dvdy(:,:,i));

        %assume that the sum of the three diaganol elements of the strain
        %rate tensor are equal to zero (by continuity)
        S33=-(S11+S22);

        %assume that the sum of the remaining terms in the strain rate tensor
        %are equal to the sum of the direcly calculated terms (5 of 9 terms
        %were calculated, explaining the '9/5' out front)
        SijSij=(9/5)*((S11.*S11)+(S12.*S12)+(S21.*S21)+(S22.*S22)+(S33.*S33));

        %compute the SGS stress using the Smagorinsky model
        epsilon_TR(:,:,i)=2*((coef*filt)^2)*(SijSij.^(3/2));

        clear S11 S12 S21 S22 S33 SijSij 
    end

    %ensemble average the time-resloved SGS dissipation
    epsilon=mean(epsilon_TR,3);
    k=k+1;
    if Data.phys.nondim
        Output.flowvarsteady(:,:,:,k)=epsilon./(charlen/(charvel^3));
        Output.varcellstdy{k}='dissipation(LES)*L/Uo^3';
    else
        Output.flowvarsteady(:,:,:,k)=epsilon;
        Output.varcellstdy{k}='dissipation(LES)';
    end
    clear term1 term2 term3 term4 term5 epsilon
    %fprintf('Done')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Complete Term Dissipation
%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.dissCOMP
    term1 = mean(Grad.dudx.^2,3);
    term2 = mean(Grad.dvdy.^2,3);
    term3 = mean(Grad.dudy.^2,3);
    term4 = mean(Grad.dvdx.^2,3);
    term5 = mean(Grad.dudy.*Grad.dvdx,3);

    if wflag == 0;
        dissipgrad = kinvisc*(2.*term1+2.*term2+3.*term3+3.*term4+2.*term5);
    else
        term6 = mean(Grad.dwdx.^2,3);
        term7 = mean(Grad.dwdy.^2,3);
        term8 = 0.5.*(term1+term2); %(dwdz)^2 & (dudz)^2 & (dudz)^2
        term9 = -0.25.*(term1+term2);
        
        dissipgrad = kinvisc*(2.*(term1 + term2              + term8) + ...
                                  term3 + term4              + 0.5.*(term1+term2) + ...
                                  term6 + 0.5.*(term1+term2) + term7 + ...
                              2.*(term5 + 2.*term9));
        clear term6 term7 term8 term9
    end
        
    k=k+1;
    if Data.phys.nondim
        Output.flowvarsteady(:,:,:,k)=dissipgrad./(charlen/(charvel^3));
        Output.varcellstdy{k}='dissipation(COMP)*L/Uo^3';
    else
        Output.flowvarsteady(:,:,:,k)=dissipgrad;
        Output.varcellstdy{k}='dissipation(COMP)';
    end
    clear term1 term2 term3 term4 term5 dissipgrad;
    
end

if Data.calc.diss1term+Data.calc.diss5term+Data.calc.dissLES + Data.calc.dissCOMP > 0
    fprintf('Done')
    T1 = etime(clock,To);
    fprintf(' %0.2i:%0.2i.%0.0f\n',floor(T1/60),floor(rem(T1,60)),...
        rem(T1,60)-floor(rem(T1,60)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vortex Identifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.calc.vortexID
    To = clock;
    fprintf('Vortex Identification...')

    if Data.calc.steady
        k = k + 1;
        cal_mat = zeros(Jdat,Idat);
        for m = 1:Jdat
            for n = 1:Idat
                velgradient=[mean(Grad.dudx(m,n),3) mean(Grad.dudy(m,n),3); ...
                    mean(Grad.dvdx(m,n),3) mean(Grad.dvdy(m,n),3)]; % velocity gradient tensor
                S=0.5*(velgradient+transpose(velgradient)); % this is symmetric
                Omega=0.5*(velgradient-transpose(velgradient)); % this is anti-symmetric

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 5. Choose computation method for vortex identification & do
                % computation
                vortex_method = {'Qcrit','Dcrit','Lamda2'};
                switch vortex_method{Data.calc.VortexIDmethod}
                    case 'Qcrit'
                        % Compares magnitude of straining & rotational motions
                        Omegan=norm(Omega,'fro'); % why normalizing these?
                        Sn=norm(S,'fro');
                        Q=0.5*(Omegan.^2-Sn.^2);
                        if Q>0
                            cal_mat(m,n)=Q;
                        end
                        Output.varcellstdy{k}='Q Critical';
                    case 'Dcrit'
                        % Cantwell's method
                        Omegan=norm(Omega,'fro');
                        Sn=norm(S,'fro');
                        Q=0.5*(Omegan.^2-Sn.^2);
                        Delta=(Q/3)^3+(det(velgradient)/2)^2;
                        if Delta>0
                            cal_mat(m,n)=Delta;
                        end
                        Output.varcellstdy{k}='D Critical';
                    case 'Lamda2'
                        % Lambda_2 method
                        T=eig(S^2+Omega^2); % principal axes?
                        Lamda_min(m,n)=min(T);%#ok
                        Lamda_max(m,n)=max(T);%#ok
                        if Lamda_max(m,n) < 0;
                            cal_mat(m,n)=Lamda_max(m,n);
                        else
                            cal_mat(m,n)=0;
                        end
                        Output.varcellstdy{k}='Lamda 2';
                end
                clear S Sn Omega Omegan Q Delta T Lamda_min Lamda_max
            end
        end
        Output.flowvarsteady(:,:,:,k)=cal_mat;
    end
    
    if Data.calc.unsteady && singlematdat
        c = c + 1;
        cal_mat = zeros(Jdat,Idat,numframes);
        for m = 1:Jdat
            for n = 1:Idat
                for p = 1:numframes
                    velgradient=[Grad.dudx(m,n,p) Grad.dudy(m,n,p); ...
                        Grad.dvdx(m,n,p) Grad.dvdy(m,n,p)]; % velocity gradient tensor
                    S=0.5*(velgradient+transpose(velgradient)); % this is symmetric
                    Omega=0.5*(velgradient-transpose(velgradient)); % this is anti-symmetric

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % 5. Choose computation method for vortex identification & do
                    % computation
                    vortex_method = {'Qcrit','Dcrit','Lamda2'};
                    switch vortex_method{Data.calc.VortexIDmethod}
                        case 'Qcrit'
                            % Compares magnitude of straining & rotational motions
                            Omegan=norm(Omega,'fro'); % why normalizing these?
                            Sn=norm(S,'fro');
                            Q=0.5*(Omegan.^2-Sn.^2);
                            if Q>0
                                cal_mat(m,n,p)=Q;
                            end
                          Output.varcellunstdy{c}='Q Critical';
                        case 'Dcrit'
                            % Cantwell's method
                            Omegan=norm(Omega,'fro');
                            Sn=norm(S,'fro');
                            Q=0.5*(Omegan.^2-Sn.^2);
                            Delta=(Q/3)^3+(det(velgradient)/2)^2;
                            if Delta>0
                                cal_mat(m,n,p)=Delta;
                            end
                            Output.varcellunstdy{c}='D Critical';
                        case 'Lamda2'
                            % Lambda_2 method
                            T=eig(S^2+Omega^2); % principal axes?
                            Lamda_min(m,n)=min(T);%#ok
                            Lamda_max(m,n)=max(T);%#ok
                            if Lamda_max(m,n) < 0;
                                cal_mat(m,n,p)=Lamda_max(m,n);
                            else
                                cal_mat(m,n,p)=0;
                            end
                            Output.varcellunstdy{c}='Lamda2';
                    end
                end
                clear S Sn Omega Omegan Q Delta T Lamda_min Lamda_max
            end
        end
          Output.flowvarunsteady(:,:,:,c)=cal_mat;
    end
    clear cal_mat

    fprintf('Done')
    T1 = etime(clock,To);
    fprintf('    %0.2i:%0.2i.%0.0f\n',floor(T1/60),floor(rem(T1,60)),...
        rem(T1,60)-floor(rem(T1,60)))
end

To = clock;
fprintf('Saving...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output Writting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.output.fileformat == 1
    if ~exist(Data.output.saveloc,'dir')
        mkdir(Data.output.saveloc)
    end
    if Data.calc.steady
        Out.flowvarsteady = Output.flowvarsteady;
        Out.varcellstdy   = Output.varcellstdy;
        %st_savename = sprintf('%s-steady-%0.0f.mat', Data.plt_info.pltname,numframes);
        st_savename = sprintf('%s-steady-%0.0f.mat', Data.output.savename,numframes);
        save(fullfile(Data.output.saveloc,st_savename),'Out','Data')
%        save([Data.output.saveloc Data.plt_info.pltname '-steady-' num2str(numframes) '.mat'],'Out','Data')
        clear Out
    end
    if Data.calc.unsteady && singlematdat
        Out.flowvarunsteady = Output.flowvarunsteady;
        Out.varcellunstdy   = Output.varcellunstdy;
        %unst_savename = sprintf('%s-unsteady-%0.0f.mat',Data.plt_info.pltname,numframes);
        unst_savename = sprintf('%s-unsteady-%0.0f.mat',Data.output.savename,numframes);
        save(fullfile(Data.output.saveloc,unst_savename),'Out','Data','-v7.3')
%        save([Data.plt_info.pltdir Data.plt_info.pltname '-unsteady-' num2str(numframes) '.mat'],'Out','Data')
        clear Out
    end
elseif Data.output.fileformat == 2
    if Data.calc.steady
        Out.flowvarsteady = Output.flowvarsteady;
        Out.varcellstdy   = Output.varcellstdy;
        %st_savename = sprintf('%s-steady-%0.0f.mat', Data.plt_info.pltname,numframes);
        st_savename = sprintf('%s-steady-%0.0f.mat', Data.output.savename,numframes);
        save(fullfile(Data.output.saveloc,st_savename),'Out','Data')
%        save([Data.plt_info.pltdir Data.plt_info.pltname '-steady-' num2str(numframes) '.mat'],'Out','Data')
        clear Out
    end
    if Data.calc.unsteady
        numvect = startframe:framestep:endframe;
        for ii = 1:numframes
            fprintf('\n')
            save_calc_multimats(Output,PLT,wflag,Data);
            Out.flowvarunsteady = permute(Output.flowvarunsteady(:,:,ii,:),[1 2 4 3]);
            Out.varcellunstdy   = Output.varcellunstdy;
            %unst_savename = sprintf('%s-unsteady-%06d.mat',Data.plt_info.pltname,numvect(ii));
            unst_savename = sprintf('%s-unsteady-%06d.mat',Data.output.savename,numvect(ii));
            save(fullfile(Data.output.saveloc,unst_savename),'Out','Data','-v7.3')
            clear Out
       end
    end
elseif Data.output.fileformat == 3
    if Data.calc.steady
        write3pifile(Output.flowvarsteady,Output.varcellstdy,time,Data.output.saveloc,[Data.plt_info.pltname,'-steady-',num2str(numframes)],Data.output.fileformat-2)
    end
    if Data.calc.unsteady && singlematdat
        write3pifile(Output.flowvarunsteady,Output.varcellunstdy,time,Data.output.saveloc,[Data.plt_info.pltname,'-unsteady-',num2str(numframes)],Data.output.fileformat-2)
    end
else
    if Data.calc.steady
        write3pifile(Output.flowvarsteady,Output.varcellstdy,time,Data.output.saveloc,[Data.plt_info.pltname,'-steady-',num2str(numframes)],Data.output.fileformat-2)
    end
    if Data.calc.unsteady
        numvect = startframe:framestep:endframe;
        write3pifile(Output.flowvarunsteady,Output.varcellunstdy,time,Data.output.saveloc,[Data.plt_info.pltname,'-unsteady-'],Data.output.fileformat-2,numvect)
    end
end
fprintf('Done')
T1 = etime(clock,To);
fprintf('                   %0.2i:%0.2i.%0.0f\n',floor(T1/60),floor(rem(T1,60)),...
    rem(T1,60)-floor(rem(T1,60)))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sub-Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,u,v,extravars,varlistnew]=read_pltmod(testname,direc,startframe,endframe,numzeros,framestep,varlist)
%
%READ_PLT Read a Tecplot plt created in FlowIQ.
%   [X,Y,U,V]=READ_PLT(TESTNAME,DIREC,STARTFRAME,ENDFRAME,NUMZEROS,FRAMESTEP)
%   reads the series of files testnameXXXX.plt from directory DIREC, where
%   XXXX is some integer beginning with STARTFRAME and ending with ENDFRAME 
%   in increments of FRAMESTEP.  XXXX has a digit length equal to NUMZEROS.  
%
%   Velocity component data are returned in U and V matrices of unknown 
%   dimension MxNxT, where M is the size in the X direction, N is the size
%   in the Y direction, and T is the total number of frames read in. 
%
%   X and Y are coordinate matrices locating the positions of the vectors
%   stored in U and V.  They have dimension MxN.
%
%   NUMZEROS and FRAMESTEP are optional, and have default values of 4 and
%   1, respectively.

%%%
%   long term plan, allow variable output argument length
%   if length(vargout)=1
%   return contents in data.(variablelist(1))=x,etc
%   maybe allow for selecting fields from input list {'X','U','Correlation'}?
%   X,Y,U,V would be default, of course?
%%%

%allow for numzeros to default to 4, and framestep to 1
if nargin<7
    varlist = {'all'};
    if nargin<6
        framestep=1;
        if nargin==4
            numzeros=4;
        end
    end
end

extravars={};


try
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numframes=endframe-startframe+1;

nameformat = sprintf('%%s%%0.%ui.dat',numzeros);
picfilename=sprintf(nameformat,testname,startframe); 

fid = fopen(fullfile(direc,picfilename));    

temp=fgetl(fid);%#ok                %TITLE="DPIV Data File"
temp=fgetl(fid);                    %VARIABLES="X" "Y" ....

%figure out what variables we have, and what their names are
[A,count,errmsg,nextindex]=sscanf(temp,'%*[^"]',1);
currindex = nextindex;
columns      = 0;
variablelist = cell(1);
var2         = cell(1);

while currindex<=length(temp)
    columns=columns+1;
    [variablelist{columns},count,errmsg,nextindex] = sscanf(temp(currindex:end),'%*["]%[^"]%*["]',1);    
    if isempty(variablelist{columns})
        currindex = currindex+1;
    else
        currindex = currindex+nextindex+1;
    end
end
r = 1;
for mm = 1:length(variablelist)
    if ~isempty(variablelist{mm})
        var2{r} = variablelist{mm};
        r = r+1;
    end
end
variablelist = var2;
columns = length(variablelist);

hh=0;varindexlist=[];varlistnew={};
if strcmp(varlist{1},'all')
    varlist=variablelist(5:end);
end
for i=1:length(varlist)
    for j=1:length(variablelist)
        if strcmp(varlist{i},variablelist{j})
            hh=hh+1;
            varindexlist(hh) = j;%#ok
            varlistnew{hh}=variablelist{j};%#ok     % in case a variable isnt typed in right
            break;
        end
    end
end

% if (nargout-4)~=length(varlist)
%     error('Number of outputs must match length of variable list')
% end

%determine size of arrays

%temp=fscanf(fid,'%*7c%u%*3c%u',2);   %Zone I=XXX J=XXX
%temp=fscanf(fid,'%*25c%u%*3c%u',2);   %Zone I=XXX J=XXX

tempstr=fgetl(fid);
 for rr=1:length(tempstr)-1
     if strcmp(tempstr(rr:rr+1),'I=')
         Iflag=rr;
     end
     if strcmp(tempstr(rr:rr+1),'J=')
         Jflag=rr;
     end
 end
Imax=sscanf(tempstr(Iflag+2:Jflag-1),'%u');
Jmax=sscanf(tempstr(Jflag+2:end),'%u');
        
fclose(fid);

%pre-allocate space for speed
u=zeros(Jmax,Imax,ceil(numframes/framestep));
v=zeros(Jmax,Imax,ceil(numframes/framestep));

numextravars = length(varindexlist);
for i = 1:numextravars
    vardata{i} = zeros(Jmax,Imax,ceil(numframes/framestep));%#ok
end

count=0;

for j=startframe:framestep:endframe
    count=count+1;

    picfilename=sprintf(nameformat,testname,j);
    fid = fopen(fullfile(direc,picfilename));
    
    temp=fgetl(fid);%#ok                    %TITLE="DPIV Data File"
    temp=fgetl(fid);%#ok                    %VARIABLES="X" "Y" ....
    temp=fgetl(fid);%#ok
%     temp=fscanf(fid,'%*7c%u%*3c%u',2);   %Zone I=XXX J=XXX
%     Imax=temp(1);
%     Jmax=temp(2);
    
    %read data for this frame
    variable=(fscanf(fid,'%g',[columns inf]))';
    numrows=size(variable,1);%#ok
    variabletemp=zeros(Jmax*Imax,columns);%#ok
%     if numrows~=Imax*Jmax               % if masking is used and not all points are printed in loaded .plt
%         tmpx=abs(diff(variable(:,1)));
%         tmpy=abs(diff(variable(:,2)));
%         indx=any(tmpx,2);
%         indy=any(tmpy,2);
%         dx=min(tmpx(indx));
%         dy=min(tmpy(indy));
%         rowindex=round((variable(:,2)-min(variable(:,2)))/dy+1);
%         colindex=round((variable(:,1)-min(variable(:,1)))/dx+1);
%         variabletemp((rowindex-1)*Jmax+colindex,:)=variable;
%         variable=variabletemp;
%     end
    
    disp(['frame ' num2str(j) ' loaded']);
    fclose(fid);

        u(:,:,count) = reshape(variable(:,3),Imax,Jmax)';
        v(:,:,count) = reshape(variable(:,4),Imax,Jmax)';
        for i=1:numextravars
            vardata{i}(:,:,count) = reshape(variable(:,varindexlist(i)),Imax,Jmax)';%#ok
        end
    
end

    x(:,:) = reshape(variable(:,1),Imax,Jmax)';
    y(:,:) = reshape(variable(:,2),Imax,Jmax)';
    for i=1:numextravars
        extravars{i} = vardata{i};%#ok
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catch err
    if fid==-1
        keyboard
        error('%s does not exist',fullfile(direc,picfilename))
    else
        rethrow(err)
    end
end

end

function [dudx,dudy,dvdx,dvdy,dwdx,dwdy] = gradient_dudxdvdy_chapra(u,v,w,dx,dy)
%
% [dudx,dudy,dvdx,dvdy] = gradient_dudxdvdy_chapra(u,v,dx,dy)
%
% calculates the vector field gradient using:
% Chapra central difference          (4th order accurate)
% Euler central difference           (2nd order accurate)
% Adams Bashforth backward Euler     (2nd order accurate)
% Adams Bashforth forward Euler      (2nd order accurate)
%
% Chris Weiland 
% Adric Eckstein
% 10.3.2007
% version 2.0

[I,J] = size(u);

dudx=zeros(I,J);
dudy=zeros(I,J);
dvdx=zeros(I,J);
dvdy=zeros(I,J);
dwdx=zeros(I,J);
dwdy=zeros(I,J);

for i=1:I
    
    %Adams Bashforth forward difference
    dudx(i,1) = (-3*u(i,1) + 4*u(i,2) - u(i,3))/(2*dx);
    dvdx(i,1) = (-3*v(i,1) + 4*v(i,2) - v(i,3))/(2*dx);
    dwdx(i,1) = (-3*w(i,1) + 4*w(i,2) - w(i,3))/(2*dx);
    
    %Adams Bashforth backward difference
    dudx(i,J) = (3*u(i,J) - 4*u(i,J-1) + u(i,J-2))/(2*dx);
    dvdx(i,J) = (3*v(i,J) - 4*v(i,J-1) + v(i,J-2))/(2*dx);
    dwdx(i,J) = (3*w(i,J) - 4*w(i,J-1) + w(i,J-2))/(2*dx);
    
    %Euler central difference
    dudx(i,2)   = (u(i,3) - u(i,1)  )/(2*dx);
    dvdx(i,2)   = (v(i,3) - v(i,1)  )/(2*dx);
    dwdx(i,2)   = (w(i,3) - w(i,1)  )/(2*dx);

    dudx(i,J-1) = (u(i,J) - u(i,J-2))/(2*dx);
    dvdx(i,J-1) = (v(i,J) - v(i,J-2))/(2*dx);
    dwdx(i,J-1) = (w(i,J) - w(i,J-2))/(2*dx);
    
    
    %Chapra central difference
    for j=3:J-2
        dudx(i,j) = -(u(i,j+2) - 8*u(i,j+1) + 8*u(i,j-1) - u(i,j-2))/(12*dx);
        dvdx(i,j) = -(v(i,j+2) - 8*v(i,j+1) + 8*v(i,j-1) - v(i,j-2))/(12*dx);
        dwdx(i,j) = -(w(i,j+2) - 8*w(i,j+1) + 8*w(i,j-1) - w(i,j-2))/(12*dx);
        
    end
    
end

for j=1:J
    
    %Adams Bashforth forward difference
    dudy(1,j) = (-3*u(1,j) + 4*u(2,j) - u(3,j))/(2*dy);
    dvdy(1,j) = (-3*v(1,j) + 4*v(2,j) - v(3,j))/(2*dy);
    dwdy(1,j) = (-3*w(1,j) + 4*w(2,j) - w(3,j))/(2*dy);
    
    %Adams Bashforth backward difference
    dudy(I,j) = (3*u(I,j) - 4*u(I-1,j) + u(I-2,j))/(2*dy);
    dvdy(I,j) = (3*v(I,j) - 4*v(I-1,j) + v(I-2,j))/(2*dy);
    dwdy(I,j) = (3*w(I,j) - 4*w(I-1,j) + w(I-2,j))/(2*dy);
    
    %Euler central difference
    dudy(2,j)   = (u(3,j) - u(1,j)  )/(2*dy);
    dvdy(I-1,j) = (v(I,j) - v(I-2,j))/(2*dy);
    dwdy(I-1,j) = (w(I,j) - w(I-2,j))/(2*dy);
    
    %Chapra central difference
    for i=3:I-2
        dudy(i,j) = -(u(i+2,j) - 8*u(i+1,j) + 8*u(i-1,j) - u(i-2,j))/(12*dy);
        dvdy(i,j) = -(v(i+2,j) - 8*v(i+1,j) + 8*v(i-1,j) - v(i-2,j))/(12*dy);
        dwdy(i,j) = -(w(i+2,j) - 8*w(i+1,j) + 8*w(i-1,j) - w(i-2,j))/(12*dy);
        
    end
    
end
end

function [dudx,dudy,dvdx,dvdy,dwdx,dwdy] = gradient_dudxdvdy_order2(u,v,w,dx,dy)
%
% [dudx,dudy,dvdx,dvdy] = gradient_dudxdvdy_order2(u,v,dx,dy)
%
% calculates the vector field gradient using:
% Euler central difference           (2nd order accurate)
% Adams Bashforth backward Euler     (2nd order accurate)
% Adams Bashforth forward Euler      (2nd order accurate)
%
% Chris Weiland 
% Adric Eckstein
% 10.3.2007
% version 2.0

[I,J] = size(u);

dudx=zeros(I,J);
dudy=zeros(I,J);
dvdx=zeros(I,J);
dvdy=zeros(I,J);
dwdx=zeros(I,J);
dwdy=zeros(I,J);

for i=1:I
    
    %Adams Bashforth forward difference
    dudx(i,1) = (-3*u(i,1) + 4*u(i,2) - u(i,3))/(2*dx);
    dvdx(i,1) = (-3*v(i,1) + 4*v(i,2) - v(i,3))/(2*dx);
    dwdx(i,1) = (-3*w(i,1) + 4*w(i,2) - w(i,3))/(2*dx);
    
    %Adams Bashforth backward difference
    dudx(i,J) = (3*u(i,J) - 4*u(i,J-1) + u(i,J-2))/(2*dx);
    dvdx(i,J) = (3*v(i,J) - 4*v(i,J-1) + v(i,J-2))/(2*dx);
    dwdx(i,J) = (3*w(i,J) - 4*w(i,J-1) + w(i,J-2))/(2*dx);
    
    %Euler central difference
    for j=2:J-1
        dudx(i,j) = (u(i,j+1) - u(i,j-1))/(2*dx);
        dvdx(i,j) = (v(i,j+1) - v(i,j-1))/(2*dx);
        dwdx(i,j) = (w(i,j+1) - w(i,j-1))/(2*dx);
        
    end
    
end

for j=1:J
    
    %Adams Bashforth forward difference
    dudy(1,j) = (-3*u(1,j) + 4*u(2,j) - u(3,j))/(2*dy);
    dvdy(1,j) = (-3*v(1,j) + 4*v(2,j) - v(3,j))/(2*dy);
    dwdy(1,j) = (-3*w(1,j) + 4*w(2,j) - w(3,j))/(2*dy);
    
    %Adams Bashforth backward difference
    dudy(I,j) = (3*u(I,j) - 4*u(I-1,j) + u(I-2,j))/(2*dy);
    dvdy(I,j) = (3*v(I,j) - 4*v(I-1,j) + v(I-2,j))/(2*dy);
    dwdy(I,j) = (3*w(I,j) - 4*w(I-1,j) + w(I-2,j))/(2*dy);
    
    %Euler central difference
    for i=2:I-1
        dudy(i,j) = (u(i+1,j) - u(i-1,j))/(2*dy);
        dvdy(i,j) = (v(i+1,j) - v(i-1,j))/(2*dy);
        dwdy(i,j) = (w(i+1,j) - w(i-1,j))/(2*dy);
        
    end
    
end
end

function [dudx,dudy,dvdx,dvdy,dwdx,dwdy] = gradient_dudxdvdy_rich4(u,v,w,dx,dy)
% 4th order richardson gradient scheme
% only does the central differences
[row,col] = size(u);

dudx=zeros(row,col);
dudy=zeros(row,col);
dvdx=zeros(row,col);
dvdy=zeros(row,col);
dwdx=zeros(row,col);
dwdy=zeros(row,col);

A=[1239 272 1036 0 -69];
k=[1 2 4 8];
for i=9:row-8
    for j=9:col-8
        cdudx=0;cdvdx=0;cdudy=0;cdvdy=0;cdwdx=0;cdwdy=0;
        for m=1:4;
            cdudx=A(m+1)*(u(i,j+k(m))-u(i,j-k(m)))/(2*k(m)*dx)+cdudx;
            cdudy=A(m+1)*(u(i+k(m),j)-u(i-k(m),j))/(2*k(m)*dy)+cdudy;
            cdvdx=A(m+1)*(v(i,j+k(m))-v(i,j-k(m)))/(2*k(m)*dx)+cdvdx;
            cdvdy=A(m+1)*(v(i+k(m),j)-v(i-k(m),j))/(2*k(m)*dy)+cdvdy;
            cdwdx=A(m+1)*(w(i+k(m),j)-w(i-k(m),j))/(2*k(m)*dx)+cdwdx;
            cdwdy=A(m+1)*(w(i+k(m),j)-w(i-k(m),j))/(2*k(m)*dy)+cdwdy;

        end
        dudx(i,j)=1/A(1)*cdudx;
        dudy(i,j)=1/A(1)*cdudy;
        dvdx(i,j)=1/A(1)*cdvdx;
        dvdy(i,j)=1/A(1)*cdvdy;
        dwdx(i,j)=1/A(1)*cdwdx;
        dwdy(i,j)=1/A(1)*cdwdy;
        
    end
end
end

function save_calc_multimats(Output,PLT,wflag,Data)

[Jdat,Idat,numf] = size(PLT.u);


for ij = 1:numf

    c = 4;
    fprintf('Unsteady Mat for Frame = %6.0f\n',ij)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity Data
    %%%%%%%%%%%%%%%%%%%%%%%%%    
    if Data.phys.nondim
        charlen = str2double(Output.phys.charlen);
        charvel = str2double(Output.phys.charvel);
        Output.flowvarunsteady(:,:,1,1) = (PLT.x)./charlen;
        Output.flowvarunsteady(:,:,1,2) = (PLT.y)./charlen;
        Output.flowvarunsteady(:,:,1,3) = (PLT.u(:,:,ij))./charvel;
        Output.flowvarunsteady(:,:,1,4) = (PLT.v(:,:,ij))./charvel;
        Output.varcellunstdy={'X/L' 'Y/L' 'U/Uo' 'V/Uo'};
    else
        Output.flowvarunsteady(:,:,1,1) = PLT.x;
        Output.flowvarunsteady(:,:,1,2) = PLT.y;
        Output.flowvarunsteady(:,:,1,3) = PLT.u(:,:,ij);
        Output.flowvarunsteady(:,:,1,4) = PLT.v(:,:,ij);
        Output.varcellunstdy={'X' 'Y' 'U' 'V'};
    end
    
    if wflag
        c = 5;
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,1,c) = PLT.w(:,:,ij)./charvel;
            Output.varcellunstdy{c}={'W/Uo'};
        else
            Output.flowvarunsteady(:,:,1,c) = PLT.w(:,:,ij);
            Output.varcellunstdy{c}={'W'};
        end
    end        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity Magnitude
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.calc.velmagvar
        if wflag == 0
            velmag=sqrt(PLT.u(:,:,ij).^2+PLT.v(:,:,ij).^2);
        else
            velmag=sqrt(PLT.u(:,:,ij).^2+PLT.v(:,:,ij).^2+PLT.w(:,:,ij).^2);
        end

        c=c+1;

        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=velmag./charvel;
            Output.varcellunstdy{c}='Velmag/Uo';
        else
            Output.flowvarunsteady(:,:,:,c)=velmag;
            Output.varcellunstdy{c}='Velmag';
        end
        clear velmag;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradients
    %%%%%%%%%%%%%%%%%%%%%%%%%
    GMETHOD = {'order2','rich4','chapra'};
    dx=PLT.x(1,2)-PLT.x(1,1);
    dy=PLT.y(2,1)-PLT.y(1,1);

    %     if ~Data.calc.gradienttype

    Grad.dudx=zeros(Jdat,Idat,'single');
    Grad.dvdx=zeros(Jdat,Idat,'single');
    Grad.dudy=zeros(Jdat,Idat,'single');
    Grad.dvdy=zeros(Jdat,Idat,'single');
    Grad.dudz=zeros(Jdat,Idat,'single');
    Grad.dvdz=zeros(Jdat,Idat,'single');
    Grad.dwdx=zeros(Jdat,Idat,'single');
    Grad.dwdy=zeros(Jdat,Idat,'single');
    
    switch lower(GMETHOD{Data.calc.gradientmethod})
        case'order2'
            if wflag > 0
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_order2(PLT.u(:,:,ij),PLT.v(:,:,ij),PLT.w(:,:,ij),dx,dy);
            else
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_order2(PLT.u(:,:,ij),PLT.v(:,:,ij),PLT.w(:,:,1),dx,dy);
            end
        case 'rich4'
            if wflag > 0
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_rich4(PLT.u(:,:,ij),PLT.v(:,:,ij),PLT.w(:,:,ij),dx,dy);
            else
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_rich4(PLT.u(:,:,ij),PLT.v(:,:,ij),PLT.w(:,:,1),dx,dy);
            end
        case 'chapra'
            if wflag > 0
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_chapra(PLT.u(:,:,ij),PLT.v(:,:,ij),PLT.w(:,:,ij),dx,dy);
            else
                [Grad.dudx(:,:),Grad.dudy(:,:),Grad.dvdx(:,:),Grad.dvdy(:,:),Grad.dwdx(:,:),Grad.dwdy(:,:)] =...
                    gradient_dudxdvdy_chapra(PLT.u(:,:,ij),PLT.v(:,:,ij),PLT.w(:,:,1),dx,dy);
            end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Vorticity
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.calc.vorticityvar
        vort=Grad.dvdx-Grad.dudy;
        c=c+1;
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=vort.*charlen/charvel;
            Output.varcellunstdy{c}='Vorticity*L/Uo';
        else
            Output.flowvarunsteady(:,:,:,c)=vort;
            Output.varcellunstdy{c}='Vorticity';
        end
        clear vort;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Strain rate, off diagonal term of 2x2 tensor
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.calc.strainratevar
        strainrate=Grad.dudy+Grad.dvdx;
        c=c+1;
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=strainrate.*charlen/charvel;
            Output.varcellunstdy{c}='Strain Rate*L/Uo';
        else
            Output.flowvarunsteady(:,:,:,c)=strainrate;
            Output.varcellunstdy{c}='Strain Rate';
        end
        clear strainrate;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradients
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.calc.gradientsvar        % gradients
        c=c+4;
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c-3)=Grad.dudx.*charlen/charvel;
            Output.flowvarunsteady(:,:,:,c-2)=Grad.dvdx.*charlen/charvel;
            Output.flowvarunsteady(:,:,:,c-1)=Grad.dudy.*charlen/charvel;
            Output.flowvarunsteady(:,:,:,c)=Grad.dvdy.*charlen/charvel;
            Output.varcellunstdy{c-3}='dudx*L/Uo';
            Output.varcellunstdy{c-2}='dvdx*L/Uo';
            Output.varcellunstdy{c-1}='dudy*L/Uo';
            Output.varcellunstdy{c}='dvdy*L/Uo';
            if wflag > 0;
                Output.flowvarunsteady(:,:,:,c+1)=Grad.dwdx.*charlen/charvel;
                Output.flowvarunsteady(:,:,:,c+2)=Grad.dwdy.*charlen/charvel;
                Output.varcellunstdy{c+1}='dwdx*L/Uo';
                Output.varcellunstdy{c+2}='dwdy*L/Uo';
                c = c+2;
            end
        else
            Output.flowvarunsteady(:,:,:,c-3)=Grad.dudx;
            Output.flowvarunsteady(:,:,:,c-2)=Grad.dvdx;
            Output.flowvarunsteady(:,:,:,c-1)=Grad.dudy;
            Output.flowvarunsteady(:,:,:,c)=Grad.dvdy;
            Output.varcellunstdy{c-3}='dudx';
            Output.varcellunstdy{c-2}='dvdx';
            Output.varcellunstdy{c-1}='dudy';
            Output.varcellunstdy{c}='dvdy';
            if wflag > 0;
                Output.flowvarunsteady(:,:,:,c+1)=Grad.dwdx;
                Output.flowvarunsteady(:,:,:,c+2)=Grad.dwdy;
                Output.varcellunstdy{c+1}='dwdx';
                Output.varcellunstdy{c+2}='dwdy';
                c = c+2;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Turbulent Kinetic Energy
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.calc.turbKEvar
        if wflag == 0;
            turbke = (PLT.u(:,:,ij)-PLT.umn).^2 + (PLT.v(:,:,ij)-PLT.vmn).^2;
        else
            turbke = (PLT.u(:,:,ij)-PLT.umn).^2 + (PLT.v(:,:,ij)-PLT.vmn).^2 + (PLT.w(:,:,ij)-PLT.wmn).^2;
        end

        c=c+1;
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=turbke./(charvel^2);
            Output.varcellunstdy{c}='TurbulentKE/Uo^2';
        else
            Output.flowvarunsteady(:,:,:,c)=turbke;
            Output.varcellunstdy{c}='TurbulentKE';
        end
        clear turbke;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Total Kinetic Energy
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.calc.totalKEvar
        if wflag == 0;
            totalke=PLT.u(:,:,ij).^2+PLT.v(:,:,ij).^2;
        else
            totalke=PLT.u(:,:,ij).^2+PLT.v(:,:,ij).^2+PLT.w(:,:,ij).^2;
        end

        c=c+1;
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c)=totalke./(charvel^2);
            Output.varcellunstdy{c}='TotalKE/Uo^2';
        else
            Output.flowvarunsteady(:,:,:,c)=totalke;
            Output.varcellunstdy{c}='TotalKE';
        end

        clear totalke;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Reynolds Stresses
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.calc.reynoldsstressvar
        c=c+3;
        if Data.phys.nondim
            Output.flowvarunsteady(:,:,:,c-2)=((PLT.u(:,:,ij)-PLT.umn).^2)./(charvel^2);
            Output.flowvarunsteady(:,:,:,c-1)=((PLT.v(:,:,ij)-PLT.vmn).^2)./(charvel^2);
            Output.flowvarunsteady(:,:,:,c)=((PLT.u(:,:,ij)-PLT.umn).*(PLT.v(:,:,ij)-PLT.vmn))./(charvel^2);
            Output.varcellunstdy{c-2}='u''^2/Uo^2';
            Output.varcellunstdy{c-1}='v''^2/Uo^2';
            Output.varcellunstdy{c}='u''v''/Uo^2';
            if wflag > 0
                Output.flowvarunsteady(:,:,:,c+1)=((PLT.w(:,:,ij)-PLT.wmn).^2)./(charvel^2);
                Output.varcellunstdy{c+1}='w''^2/Uo^2';
                c=c+1;
            end
        else
            Output.flowvarunsteady(:,:,:,c-2)=(PLT.u(:,:,ij)-PLT.umn).^2;
            Output.flowvarunsteady(:,:,:,c-1)=(PLT.v(:,:,ij)-PLT.vmn).^2;
            Output.flowvarunsteady(:,:,:,c)=(PLT.u(:,:,ij)-PLT.umn).*(PLT.v(:,:,ij)-PLT.vmn);
            Output.varcellunstdy{c-2}='u''^2';
            Output.varcellunstdy{c-1}='v''^2';
            Output.varcellunstdy{c}='u''v''';
            if wflag > 0
                Output.flowvarunsteady(:,:,:,c+1)=(PLT.w(:,:,ij)-PLT.wmn).^2;
                Output.varcellunstdy{c+1}='w''^2';
                c=c+1;
            end

        end
        clear reystressuprime reystressvprime reystresswprime reystressupvp reystressupwp reystressvpwp;
    end

    Out.flowvarunsteady = Output.flowvarunsteady;
    Out.varcellunstdy   = Output.varcellunstdy;
    unst_savename = sprintf('%s%s-unsteady-%06d.mat',Data.output.saveloc, Data.plt_info.pltname,ij);
    save(unst_savename,'Out','Data')
    clear Out

end
end

function write3pifile(flowvarunsteady,varcellunstdy,time,direcplt,nameplt,fileformat,framenum)
% this function writes an unsteady 3pi file.  the calculations are in
% flowvarunsteady (Jdat x Idat x numberofframes x numberofvariables).  the files
% are stored in direcplt.  fileformat=1(.dat file) or =2 (.mat file).
% flowvariablelist is a logical array (1=variable active, 0=variable
% inactive).  time is a vector containing the time at each frame (length
% numframes).  x and y are matrices.  u and v are 


[Jdat,Idat,numframes,numvars]=size(flowvarunsteady);

varstring=[];
for j=1:length(varcellunstdy)     % append cell array into single character array
    varstring=[varstring '"' varcellunstdy{j} '", '];%#ok
end
varstring=varstring(1:end-2);       % take off final comma and space

if fileformat==1 %|| fileformat==3      % tecplot .dat formatting, block
    numstring=deblank(repmat('%+11.3e\t',1,Idat));     % numstring should have a total of Idat formatting entries
    
    filename=[nameplt '.dat'];

    fid = fopen([direcplt filename], 'wt');

    fprintf(fid,'TITLE = "DPIV Data File"\n');            % Block formatted tecplot .dat file
    fprintf(fid,['VARIABLES =' varstring '\n']);
    tempmat=zeros(Jdat,Idat*numvars)';%#ok

    for k=1:numframes
        fprintf(fid,['Zone T= "' num2str(time(k)) '", I=%3.0f J=%3.0f F=BLOCK C=BLACK SOLUTIONTIME=' num2str(time(k)) '\n'],Idat,Jdat);      % write new zone text

        for j=1:numvars
            tempmat=flowvarunsteady(:,:,k,j)';
            fprintf(fid,[numstring '\n'],tempmat);
        end
    end
    fclose(fid);

elseif fileformat==2
    
    for k=1:numframes
        filename = [nameplt,sprintf('%06.0f.dat',framenum(k))];
        fid = fopen([direcplt filename], 'wt');
        
        fprintf(fid,'TITLE = "DPIV Data File"\n');            % Point formatted tecplot .dat file
        fprintf(fid,['VARIABLES =' varstring '\n']);
        tempmat=zeros(Jdat,Idat*numvars)';%#ok
        
        fprintf(fid,['Zone T= "' num2str(time(k)) '", I=%3.0f J=%3.0f F=POINT C=BLACK SOLUTIONTIME=' num2str(time(k)) '\n'],Idat,Jdat);      % write new zone text
        for j=1:Jdat
            for i=1:Idat
                for x=1:numvars
                    fprintf(fid,'%14.7g',flowvarunsteady(j,i,k,x));
                end
                fprintf(fid,'\n');
            end
        end
        
        fclose(fid);
    end
else        % matlab output format
    save([direcplt nameplt],'varcellunstdy','flowvarunsteady');
%     save(nameplt,'varcellunstdy','flowvarunsteady','-v7.3');%Outputs into the 3pi folder
end
end