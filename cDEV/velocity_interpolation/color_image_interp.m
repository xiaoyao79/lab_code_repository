function [data] = color_image_interp(data,info)
% Subfunction for (1) extracting color image information from the data
% and (2) interpolating from color data to velocity values

% Cut out ROI from Color Images
data.imcolor = zeros([info.s1 info.s2 info.s3 info.s4]);

for ii = 1:info.s4
    data.imcolor(:,:,:,ii) = double(data.imcdopp(:,:,:,ii)).*...
        double(repmat(data.BWreg,[1 1 info.s3]));
end

% Subtract Grayscale from Image
Iblank = zeros(1,1,1);
data.imgray = zeros([info.s1 info.s2 info.s3 info.s4]);

for ii = 1:info.s4
    Icolor = data.imcolor(:,:,:,ii);
    se1 = strel('disk',5);
    se2 = strel('disk',7);
    
    % Identify Red-Blue
    temp = abs(Icolor(:,:,1)-Icolor(:,:,3));
    temp = temp(temp~=0);
    
    Itemp  = ones(size(Icolor(:,:,1)));
    Itemp(abs(Icolor(:,:,1)-Icolor(:,:,3)) <= (median(temp(:))+std(temp(:)))) = Iblank;
    Itemp = imdilate(imerode(Itemp,se1),se2);
    slug1 = Icolor.*repmat(Itemp,[1 1 3]);
    
    % Identify Red-Green
    temp = abs(Icolor(:,:,1)-Icolor(:,:,2));
    temp = temp(temp~=0);
    
    Itemp  = ones(size(Icolor(:,:,1)));
    Itemp(abs(Icolor(:,:,1)-Icolor(:,:,2)) <= (median(temp(:))+std(temp(:)))) = Iblank;
    Itemp = imdilate(imerode(Itemp,se1),se1);
    slug2 = Icolor.*repmat(Itemp,[1 1 3]);
    
    % Identify Green-Blue
    temp = abs(Icolor(:,:,2)-Icolor(:,:,3));
    temp = temp(temp~=0);
    
    Itemp  = ones(size(Icolor(:,:,1)));
    Itemp(abs(Icolor(:,:,2)-Icolor(:,:,3)) <= (median(temp(:))+std(temp(:)))) = Iblank;
    Itemp = imdilate(imerode(Itemp,se1),se1);
    slug3 = Icolor.*repmat(Itemp,[1 1 3]);
    
    slug1(slug1 ~= 0) = 1;
    slug2(slug2 ~= 0) = 1;
    slug3(slug3 ~= 0) = 1;
    
    slug = (slug1 + slug2 + slug3);
    slug(slug ~= 0) = 1;
    Icolor = Icolor.*slug;
    
    figure(20); imshow(uint8(Icolor));
    set(gcf,'Position',[100 100 0.6*info.s2 0.6*info.s1]);
    pause(1e-1);
    
    data.imgray(:,:,:,ii)  = data.imcolor(:,:,:,ii) - Icolor;
    data.imcolor(:,:,:,ii) = Icolor;
end

% Select ColorScale for De-aliasing
figure(1);
set(gcf,'Position',[100 100 0.6*info.s2 0.6*info.s1]);
imagesc(uint8(data.imcdopp(:,:,:,1)));

Vmin = input('Enter a Value for Vmin \n');
Vmax = input('Enter a Value for Vmax \n');

% Crop Color Scale From Images for Dealiasing
[~,~,I_cscale,~]=imcrop(uint8(data.imcdopp(:,:,:,1)));  
fprintf('Select Color Scale Region for Velocity Interpolation \r')
set(gcf,'Position',[100 100 0.6*info.s2 0.6*info.s1]);

image(I_cscale);
set(gcf,'Position',[100 100 0.6*info.s2 0.6*info.s1]);
scale3 = permute(mean(I_cscale, 2),[1,3,2]);
full_scale = linspace(Vmax,Vmin,size(scale3,1))';

for ii = 1:info.s4
    [height, width,~] = size(...
        data.imcolor(data.ystart:data.yend,data.xstart:data.xend,:,ii));
    Iflat = reshape(...
        data.imcolor(data.ystart:data.yend,data.xstart:data.xend,:,ii),...
        [height*width,3]);
    
    figure(4);
    [~, min_index] = min(abs(full_scale));
    
    %spline fit for colorscale
    knots = augknt([Vmin,Vmin*3/4,Vmin/2,Vmin*1/4,0,Vmax*1/4,Vmax/2,Vmax*3/4,Vmax],4,[1,1,1,3,1,1,1]);
    ww    = ones(size(full_scale.'));
    ww(1) = 10; ww(end) = 10; ww(min_index) = 10;
    sp_rgb = spap2(knots,4,full_scale.',[scale3(:,1),scale3(:,2),scale3(:,3)]',ww);
    sp_rgb.coefs;
    fnplt(sp_rgb)
    
    hold on
    plot3(Iflat(1:10:end,1),Iflat(1:10:end,2),Iflat(1:10:end,3),'ro');
    plot3(scale3(:,1),scale3(:,2),scale3(:,3),'k');
    hold off
    xlabel('red'),ylabel('green'),zlabel('blue')
    legend('reconstructed color scale','original sampled scale','measured values')
    rgb_error = inline('norm(fnval(sp,x)-rgb)','x','sp','rgb');
    
    pause(1E-1);
    % Run interpolation
%     [data.UU(:,:,ii),data.velocity_map(:,:,ii),...
%         data.reconstructed(:,:,:,ii),data.mapp] = ...
%         interpolation(...
%         data.imcolor(data.ystart:data.yend,data.xstart:data.xend,:,ii),...
%         rgb_error, Vmin, Vmax,sp_rgb,full_scale);
    
%%
    [data.UU(:,:,ii),data.velocity_map(:,:,ii),...
        data.reconstructed(:,:,:,ii),data.mapp]  = ...
    colorreg_interp(data.imcolor(data.ystart:data.yend,data.xstart:data.xend,:,ii),...
        sp_rgb,full_scale,I_cscale);
end
data.Vmin = Vmin;
data.Vmax = Vmax;