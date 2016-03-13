function [data] = window_data(data,info)
% Subfunction for generating bulk window and subregion windows, then apply
% subregion window to limit the amount of necessary data for evaluation

% Block for generating bulk window mask (taken from prana)
data.mask   = ones([info.s1 info.s2]);
roiwindow   = CROIEditor((uint8(data.imcdopp(:,:,:,1))));

while isempty(roiwindow.labels)
    addlistener(roiwindow,'MaskDefined',@your_roi_defined_callback);
    drawnow
end

data.mask(roiwindow.labels==0) = 0;
close('Analyzer - ROI Editor')

% Sample bounded window
imbound = zeros([info.s1 info.s2 info.s3 info.s4]);
for ii = 1:info.s4
    for i = 1:info.s1
        for j = 1:info.s2
            if data.imcdopp(i,j,1,ii) >= 0 &&...
                    data.imcdopp(i,j,2,ii) >= 0 && ...
                    data.imcdopp(i,j,3,ii) >= 0
                imbound(i,j,:,ii) = 1;
            end
        end
    end
end

imbound = median(mean(imbound,4),3).*data.mask;
imbound(imbound < 1) = 0;

% Generate subregion window mask
se = strel('disk',2);
% se1= strel('rectangle',[9 11]);
se1= strel('disk',5);

BWReg = imfill(imclose(imdilate(imbound,se),se1));

% Generate subregion window mask
se = strel('disk',2);
% se1= strel('rectangle',[9 11]);
se1= strel('disk',2);

BWReg = imdilate(imerode(BWReg,se),se1);
data.BWreg = double(bwselect(BWReg,data.xloc,data.yloc));

% Find region start and end locations
[r,c] = find(data.BWreg == 1);
data.xstart = min(c(:))-1;
data.xend   = max(c(:))+1;
data.ystart = min(r(:))-1;
data.yend   = max(r(:))+1;

% Window Theta Matrix & Rho Matrix
data.Theta = data.Theta.*data.BWreg;
data.Rho   = data.Rho.*data.BWreg;