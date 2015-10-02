scale = round(0.3*max(im(:)));
im1   = im; im1(im1 < scale) = 0;

imhess = hessian2(double(im1).*double(mask),[3 9]);

im_mask1 = im2bw(uint8(double(intmax('uint8'))/max(imhess(:))*imhess),0.2);
se = strel('disk',15);
im_mask2 = imclose(im_mask1,se);
im_mask3 = bwmorph(im_mask2,'shrink',50);
%%
for i = 2:size(im_mask3,2)-1
    for j = 2:size(im_mask3,1)-1
        if im_mask3(j,i) == 1
            const = im_mask3(j,i);
            tmat  = im_mask3(j-1:j+1,i-1:i+1);
            mmat  = reshape(tmat,[size(tmat,2)*size(tmat,1) 1]);
            eval  = sum(mmat)-const;
            if eval == 0
                im_mask3(j,i) = 0;
            end
        end 
    end
end
%%
im_mask4 = imdilate(im_mask3,se);