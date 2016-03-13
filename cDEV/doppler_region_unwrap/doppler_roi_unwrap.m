function [data] = doppler_roi_unwrap(data,info)
% Set Up Grids and Unwrap ROIs
temp1 = zeros([numel(data.yloc) 1]); 
temp2 = temp1;

if max(data.Theta(:)) == 0
    data.Theta = data.THETA;
end
if max(data.Rho(:)) == 0
    data.Rho = data.RHO;
end

for ii = 1:numel(data.yloc)
    temp1(ii) = data.Theta(round(data.yloc(ii,1)),round(data.xloc(ii,1)));
    temp2(ii) = data.Rho(round(data.yloc(ii,1)),round(data.xloc(ii,1)));
end

% Select two tick marks seperated by 1cm
temp = figure(1); imagesc(uint8(data.imcdopp(:,:,:,1)));
set(gcf,'position',[100 100 0.6*info.s2 0.6*info.s1]);
[data.drx,data.dry] = getpts(temp);

data.dR = 1/sqrt((data.drx(2)-data.drx(1))^2+(data.dry(2)-data.dry(1))^2);

minT = min(temp1)-10;
maxT = max(temp1)+10;

% minR = min(temp2);
maxR = max(temp2);

% minT = min(data.Theta(:));
% maxT = max(data.Theta(:));

minR = min(data.Rho(data.Rho~=0))+20;

% maxR = max(data.Rho(:));

% ind1 = find(temp1 == minT);
% ind2 = find(temp1 == maxT);
% ind3 = find(temp2 == minR);
% ind4 = find(temp2 == maxR);

Xtform = repmat(linspace(deg2rad(maxT),deg2rad(minT),size(data.UU,2)),[size(data.UU,1) 1]);
Ytform = repmat(linspace(minR,maxR,size(data.UU,1))',[1 size(data.UU,2)])*data.dR;

data.T = data.THETA(data.ystart:data.yend,data.xstart:data.xend);
data.R = data.RHO(data.ystart:data.yend,data.xstart:data.xend)*data.dR;

for ii = 1:info.s4
    data.Vtform(:,:,ii) = griddata(data.T,data.R,data.changemap(:,:,ii),...
        Xtform,Ytform,'natural');
    data.Gtform(:,:,ii) = griddata(data.T,data.R,...
        squeeze(mean(data.imgray(data.ystart:data.yend,data.xstart:data.xend,:,ii),3)),...
        Xtform,Ytform,'natural');
end
data.Xtform = Xtform;
data.Ytform = Ytform;