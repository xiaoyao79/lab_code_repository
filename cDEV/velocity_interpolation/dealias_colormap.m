function [changemap] = dealias_colormap(velocity_map, Vmin, Vmax, mapp)
% CMM Version
% function [change_map] = dealias_colormap(velocity_map, tt, xx, Vmin, Vmax, mapp, x_m, first_t,last_t, t_E, x_m_loc)
%
% 2D Doppler
% function [change_map] = dealias_colormap(velocity_map, Vmin, Vmax, mapp)
%
% Note: This is not the proper format for commenting code
%
% Version 1.0
% Written by Brett Meyers on 08/17/2014, officially published at 7:45 PM

% Initialize certain properties
vmap = velocity_map;
vmap(abs(vmap) < min([abs(mean(velocity_map(velocity_map>0))) abs(mean(velocity_map(velocity_map<0)))])) = 0;
% vmap(abs(vmap) < max([abs(median(velocity_map(velocity_map>0))) abs(median(velocity_map(velocity_map<0)))])) = 0;

% vmap(abs(vmap) < median(abs(vmap(vmap~=0)))) = 0;
% keyboard
% vmap(abs(vmap) < max([median(abs(vmap(vmap>0))) median(abs(vmap(vmap<0)))])) = 0;
% keyboard
binary_mask = velocity_map;         % This initializes the binary mask
binary_mask(binary_mask ~=0) = 1;   % This binary mask is generated for any needed filtering of spurious data

bw_msk = sign(vmap);        % This seperates the flow field in to negative, positive, and zero signage
bw_msk(bw_msk ==  1) = 0;           % Positive data is not needed for labeling, in positive data should not be aliased
bw_msk(bw_msk == -1) = 1;           % This generates black/white binary regions where negative velocites exist

% keyboard
L = bwlabel(bw_msk,4);              % Label regions where negative regions exist w/ connectivity
L = L;                              % This bumps up the levels to identify all zeros and positive data as a region
L = L.*bw_msk;                      % Binary filtering is used to ensure zero data is not included in the positive level

bw_msk1 = sign(vmap);        % This seperates the flow field in to negative, positive, and zero signage
bw_msk1(bw_msk1 == -1) = 0;           % Positive data is not needed for labeling, in positive data should not be aliased
bw_msk1(bw_msk1 ==  1) = 1;           % This generates black/white binary regions where negative velocites exist

L1 = bwlabel(bw_msk1,4);             % Label regions where negative regions exist w/ connectivity
L1 = L1+max(L(:));                  % This bumps up the levels to identify all zeros and positive data as a region
L1 = L1.*bw_msk1;
L = L+L1;

med_map = zeros(size(velocity_map));    % This initializes the median value color scan

% This for-loop goes through each connected region and calculates the
% median value of all measurements
for t = 1:max(L(:))
    [r,c]   = find(L == t);
    for s = 1:size(r,1)
        val_stack(s,1) = velocity_map(r(s),c(s));
    end
    med_val(t,1) = median(val_stack);
    clear val_stack
end

% This section creates the colorized median value scan
M = L;                          % This initializes the median value map for reference value replacement
for t = 1:max(L(:))
    M(M==t) = med_val(t,1);     % This replaces region values with median values
end

figure(1000); imagesc(M); colormap(mapp); caxis([Vmin Vmax]);

% This section determines which regions are connected to one-another
Mbw = M;                    % Initialize black/white binary map for region pixel indexing
Mbw(Mbw~=0) = 1;            % Binarization by setting all regions not equal to zero to 1
CC = bwconncomp(Mbw,4);     % Identify connected regions in the image, returns structure with pixel index (CC.PixelIndxList)

% This section will go through and determine where there are aliased
% measurements in each region, and replace through de-aliasing
vmap_dealias = zeros(size(velocity_map));
for t = 1:size(CC.PixelIdxList,2)
    Zbw = zeros(size(velocity_map));                    % This initializes a zero black/white binary map
    Zbw(CC.PixelIdxList{t}) = L(CC.PixelIdxList{t});    % Region-identified pixel indicies for current grouping are placed in map
    
    % If there is only one region in the group, there is no need to deliase
    if numel(unique(Zbw(Zbw>0))) == 1
        Zbw1 = zeros(size(Zbw));
        Zbw1(Zbw > 0) = 1;
        vmap_dealias = vmap_dealias+velocity_map.*Zbw1;  % Add region points to list
        
        % This is where multiple regions are present, the loops identify aliased regions and attempt to de-alias
    else
        lst = unique(Zbw(Zbw>0));   % List the number of unique regions in the map
        
        % Determine the number of points wihtin each region
        for s = 1:size(lst,1)
            cnt(s,1) = numel(Zbw(Zbw == lst(s,1)));
        end
        
        % Theory dictates that the region with the greatest number of points is
        % the region with the proper set of velocity measurements.  Other
        % regions in the flow with a different sign are aliased.  De-alias
        % these regions.
        %         keyboard
        gate = find(cnt == max(cnt(:)));
        sgn  = sign(med_val(lst(gate),1));
        
        % This is where deliasing is done
        if sgn == 1
            for s = 1:size(lst,1)
                if s == gate
                    Zbw1 = zeros(size(Zbw));    % Initialize a place-holder mask
                    Zbw1(Zbw == lst(s,1)) = 1;  % Binarize kept points to 1
                    
                    % Add preserved points (at base value) to de-aliased map
                    vmap_dealias = vmap_dealias+velocity_map.*Zbw1;
                else
                    if sign(med_val(lst(s))) == sgn
                        Zbw1 = zeros(size(Zbw));
                        Zbw1(Zbw == lst(s,1)) = 1;
                        % Add other preserved points (at base value) to de-aliased map
                        vmap_dealias = vmap_dealias+velocity_map.*Zbw1;
                    else
                        % This is the section that will perform all dealiasing procedures
                        [r c] = find(Zbw == lst(s,1));      % Find the index values for every point in the aliased region
                        vmp = zeros(size(velocity_map));    % Initialize a place-holder map
                        
                        % De-alias velocity measurement values
                        %                         if abs(med_val(lst(s))) > med_val(lst(gate))
                        for q = 1:size(r,1)
                            %                             if abs(velocity_map(r(q),c(q))) > med_val(lst(gate))
                            %                                 keyboard
                            vmp(r(q),c(q)) = velocity_map(r(q),c(q))+(Vmax-Vmin);
                            %                             end
                        end
                        %                         end
                        % Add de-aliased values to de-aliased map
                        vmap_dealias = vmap_dealias+vmp;
                    end
                end
            end
        elseif sgn == -1
            for s = 1:size(lst,1)
                if s == gate
                    Zbw1 = zeros(size(Zbw));
                    Zbw1(Zbw == lst(s,1)) = 1;
                    vmap_dealias = vmap_dealias+velocity_map.*Zbw1;
                else
                    if sign(med_val(lst(s))) == sgn
                        Zbw1 = zeros(size(Zbw));
                        Zbw1(Zbw == lst(s,1)) = 1;
                        vmap_dealias = vmap_dealias+velocity_map.*Zbw1;
                    else
                        % This is the section that will perform all dealiasing
                        % procedures
                        [r c] = find(Zbw == lst(s,1));
                        vmp = zeros(size(velocity_map));
                        
                        %                     if med_val(lst(s,1)) > abs(med_val(lst(gate)))
                        for q = 1:size(r,1)
                            %                     if (velocity_map(r(q),c(q))) > abs(med_val(lst(gate)))
                            vmp(r(q),c(q)) = velocity_map(r(q),c(q))+(Vmin-Vmax);
                            %                     end
                        end
                        %                     end
                        vmap_dealias = vmap_dealias+vmp;
                    end
                end
            end
        end
        clear cnt
    end
end

for i = 1:size(velocity_map,2)
    for j = 1:size(velocity_map,1)
        if vmap_dealias(j,i) == 0;
            vmap_dealias(j,i) = velocity_map(j,i);
        end
    end
end

for i = 1:size(velocity_map,2)
    for j = 1:size(velocity_map,1)
        if vmap_dealias(j,i) > 0.9*(2*Vmax) || vmap_dealias(j,i) < 0.9*(2*Vmin);
            vmap_dealias(j,i) = velocity_map(j,i);
        end
    end
end

% keyboard
%
vmp_smooth = vmap_dealias;
%
% keyboard
% Relics from old code, unsure of their importance but just trying to
% ensure everything will run properly
% t=tt;
% x = flipdim(xx,2);
%
% dx = (x(2)-x(1))/100;               %% Define dx
% dt = (t(2)-t(1));
%
% [T,X] = meshgrid(t,x);

velocity_map2 = vmp_smooth;
velocity_map3 = velocity_map2;

% Check first_t and correct if within the wrong range
% if first_t<4  || first_t>max(t_E) || last_t<=first_t
%     V_mv=mean(velocity_map2(round(x_m_loc)-35:round(x_m_loc)-25,:));
%     V_mv=smooth(V_mv,15);
%
%     %Find first and last t
%     y=V_mv;
%     [ind,peaks] = findpeaks(-V_mv);
%     t_E_first=min(t_E);
%     first_t=max(ind(find((t_E_first-ind)>0)));
%     last_t=min(ind(find((t_E_first-ind)<0)));
% end
%
% if first_t==0;
%     first_t=1;
% end

% Create Umap which is velocity at mitral valve inlet over E wave
% Umap=zeros(size(velocity_map2,1),size(velocity_map2,2));
% x_mitral=size(velocity_map,1)-x_m;
%
% x_upper=round(x_mitral)+(7);
% x_lower=round(x_mitral)-(35);
%
% column=[first_t last_t  last_t first_t];
% row=[x_upper x_upper x_lower x_lower];
% Iroi2=roipoly(size(velocity_map2,1),size(velocity_map2,2),column,row);
%
% Umap = Umap + Iroi2.*(velocity_map2);
%
% % Determine Um value for isovelocity contours
% Umap_block=Umap(round(x_m_loc)-25:round(x_m_loc)-10,first_t:last_t)';
%
% Umap_med=max(Umap_block);
% figure(2),plot(Umap_med)
% c=1;
% for i=1:length(Umap_med)
%     if Umap_med(i)<mean(Umap_med)+1.5*std(Umap_med)&&Umap_med(i)>mean(Umap_med)-1.5*std(Umap_med)
%         Umap_med2(c)=Umap_med(i);
%         c=c+1;
%     end
% end
%
% Umap_med_smooth=smooth(Umap_med2,3);
% max(Umap_med_smooth)
% Um=max(Umap_med_smooth);
% keyboard
% vmap_dealias1 = imfill(im2bw(vmap_dealias>min([median(abs(vmap(vmap>0))) median(abs(vmap(vmap<0)))])),'holes').*velocity_map2;

vmap_dealias1 = imfill(im2bw(vmap_dealias > mean(abs(vmap(vmap~=0)))),'holes').*velocity_map2;
% vmap_dealias1 = imfill(im2bw(vmap_dealias>0),'holes').*velocity_map2;
for k=1:size(velocity_map,2);
    for i=1:size(velocity_map,1);
        if vmap_dealias1(i,k)~=0;
            if vmap_dealias1(i,k) < 0
                vmap_dealias1(i,k)=vmap_dealias1(i,k)+(Vmax-Vmin);
            end
        end
    end
end

% vmap_dealias2 = imfill(im2bw(-vmap_dealias>min([median(abs(vmap(vmap>0))) median(abs(vmap(vmap<0)))])),'holes').*velocity_map2;
vmap_dealias2 = imfill(im2bw(-vmap_dealias > mean(abs(vmap(vmap~=0)))),'holes').*velocity_map2;
% vmap_dealias2 = imfill(im2bw(-vmap_dealias>0),'holes').*velocity_map2;
for k=1:size(velocity_map,2);
    for i=1:size(velocity_map,1);
        if vmap_dealias(i,k)~=0;
            if vmap_dealias2(i,k) > 0
                vmap_dealias2(i,k)=vmap_dealias2(i,k)+(Vmin-Vmax);
            end
        end
    end
end

vmap_dealias3 = vmap_dealias1 + vmap_dealias2;

for i = 1:size(velocity_map,2)
    for j = 1:size(velocity_map,1)
        if vmap_dealias3(j,i) == 0;
            vmap_dealias3(j,i) = velocity_map(j,i);
        end
    end
end

% keyboard
velocity_map2 = vmap_dealias;
% velocity_map2 = smoothn(vmap_dealias,'robust');
% velocity_map2 = smoothn(velocity_map2,20);
changemap = velocity_map2;