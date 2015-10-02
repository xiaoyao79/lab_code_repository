function [VCP,VCN,VLP,VLN,CircAP,CircAN,CircLP,CircLN] = vortex_circ_main(U,V,cal_mat,mask,vorticity,dx,swirl_thresh,a_filt,vector_spacing,c_meth)
%%
% swirl_tresh -  percent (0 to 1) of max to draw vortex borders
% a_filt      -  ignore vortices smaller than a specific area
% vector spacing - spacing on grid between points
% c_meth      - Contour method for identifying vortices, either WeightC or
% WeightW (Area by Cal_mat vs Area by Vorticity)
%
% NV is the number of vortices identified
% VC*         - [NV,2], stored as (i,j) pairs
% CircL*      - [1,NV], circulation in cm^2/s (?) using line integral around perimeter
% VA*         - [1,NV], area of each vortex in cm^2 (?)
% cal_mat_new - modified cal_mat swirl strength, not sure what the point is
% CircA*      - [1,NV], circulation in cm^2/s (?) using area integral of vorticity
%
% percent of max to draw vortex borders
% swirl_thresh = 0.025;
% a_filt = 0.0005;    %ignore vortices smaller than a_filt (in cm^2)
% vector_spacing = 1; %because vectors are sampled at 1 vector/pix
% c_meth ='weightC';

% MASK = repmat(mask,[1 1 size(cal_mat,3)]);
MASK = mask;

tic
fprintf('Frame %i : ',size(cal_mat,3))

%Find positive vortices
% VCP - Vortex Centers, Positive
% CircLP - Positive Line Circulation
% VAP - Vortex Areas, Positive
% cal_mat_new - IDK
% CircAP - Positive Area Circulation

[VCP, CircLP, VAP, ~, CircAP] = ...
    vortex_id(cal_mat', swirl_thresh, U',V',...
    dx, a_filt, vector_spacing, vorticity',c_meth);

VCP = round(VCP);

%empty matrix, 1 marks that vortex as valid, 0 not valid
good_vortex_list = zeros(size(CircAP));

%modified from kelley's MRI_vortex_id_CS.m
for n=1:size(VCP,1)
    %if the location is non-zero?
    if VCP(n,1)>0;
        %is the center within the heartmask for LV?
        if MASK(VCP(n,1),VCP(n,2))==1
            %mark current vortex valid in list
            good_vortex_list(n) = 1;
        end
    end
end

%figure out which valid vortex is most positive (CW)
[strength, n] = max(CircAP.*good_vortex_list);

%no valid positive vortex this timestep
if strength<=0
    VLP = [1 1 0 0].';
else
    VLP = [VCP(n,1), VCP(n,2), CircAP(1,n), VAP(1,n) ].';
end

%Find negative vortices
% VCN - Vortex Centers, Positive
% CircLN - Positive Line Circulation
% VAN - Vortex Areas, Positive
% cal_mat_new - IDK
% CircAN - Positive Area Circulation
[VCN, CircLN, VAN, ~, CircAN] = ...
    vortex_id(cal_mat.', -swirl_thresh, U.',...
    V.',dx, a_filt, vector_spacing, vorticity.',c_meth);

VCN = round(VCN);

%empty matrix, 1 marks that vortex as valid, 0 not valid
good_vortex_list = zeros(size(CircAN));

%modified from kelley's MRI_vortex_id_CS.m
for n=1:size(VCN,1)
    %if the location is non-zero?
    if VCN(n,1)>0;
        %is the center within the heartmask for LV?
        if mask(VCN(n,1),VCN(n,2))==1
            %mark current vortex valid in list
            good_vortex_list(n) = 1;
        end
    end
end

%figure out which valid vortex is most negative (CCW)
[strength, n] = min(CircAN.*good_vortex_list);
if strength>=0  %no valid negative vortex this timestep
    VLN = [1 1 0 0].';
else
    VLN = [VCN(n,1), VCN(n,2), CircAN(1,n), VAN(1,n) ].';
end
keyboard
fprintf('%.3fs\n',toc)