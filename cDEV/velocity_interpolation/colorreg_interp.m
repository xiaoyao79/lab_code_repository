function [UU,velocity_map,reconstructed,mapp] = ...
    colorreg_interp(imgcolor,sp_rgb,full_scale,I_cscale)

velocity_map = zeros([size(imgcolor,1)*size(imgcolor,2) 1]);
colorimg = reshape(imgcolor,[size(imgcolor,1)*size(imgcolor,2) size(imgcolor,3)]);

tic
for i = 1:size(colorimg,1)
    if mean(colorimg(i,:)) ~= 0
%         diffcolor = repmat(squeeze(imgcolor(i,j,:))',[size(I_cscale,1) 1]) -...
%                 squeeze(double(median(I_cscale,2)));
        try    
        diffcolor = repmat(colorimg(i,:), [size(I_cscale,1) 1]) -... 
            squeeze(double(median(I_cscale,2)));
        cdiffvect = sum(abs(diffcolor),2);
        
        mincdiff = find(cdiffvect == min(cdiffvect(:)));
        
        
            velocity_map(i,:) = full_scale(mincdiff(1));
        catch
            keyboard
        end
    else
        
    end
end
velocity_map = reshape(velocity_map,[size(imgcolor,1) size(imgcolor,2)]);
toc

% for i=1:size(imgcolor,1)
%     for j=1:size(imgcolor,2)
%         if mean(imgcolor(i,j,:)) ~= 0
%             diffcolor = repmat(squeeze(imgcolor(i,j,:))',[size(I_cscale,1) 1]) -...
%                 squeeze(double(median(I_cscale,2)));
%             
%             cdiffvect = sum(abs(diffcolor),2);
%             
%             mincdiff = find(cdiffvect == min(cdiffvect(:)));
%             
%             try
%                 velocity_map(i,j) = full_scale(mincdiff(1));
%             catch
%                 keyboard
%             end
%         else
%             
%         end
%     end
% end
% toc

% go back from velocities to colors
reconstructed = permute(uint8(fnval(sp_rgb,velocity_map)),[2,3,1]);

%plot some results
 mapp = flipud(fnval(sp_rgb,full_scale')')./255;
for i=1:length(full_scale)
    for j=1:3,
        if mapp(i,j) > 1
            mapp(i,j) = 1;
        end,
        if mapp(i,j) < 0
            mapp(i,j) = 0;
        end
    end
end

for i=1:length(full_scale)
    for j=1:3
        if mapp(i,j) > 1
            mapp(i,j) = 1;
        end
        if mapp(i,j) < 0
            mapp(i,j) = 0;
        end
    end
end

UU=flipdim(velocity_map,1);