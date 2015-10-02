% Set current directory to Mouse Data Analysis M-files Folder
clear maindir
maindir = ('/Volumes/Freki/AEThER_Lab_research/AAA_data/analysis_mfiles');

% Add seperated folder directories for additional M-files needed in
% processing
addpath(fullfile(maindir,'color_scale_interp'));
addpath(fullfile(maindir,'dealias_colormap'));
addpath(fullfile(maindir,'smoothing'));
addpath(fullfile(maindir,'discrete_difference'));
addpath(fullfile(maindir,'outlier_detection'));
addpath(fullfile(maindir,'POD'));
addpath(fullfile(maindir,'Validation'));
%%
x = data.x;
y = data.y;
u = data.U;
v = repmat(data.v,[1 1 size(data.U,3)]);

xp = x/data.channelhgt;
yp = y/data.channelhgt;
up = u/data.umean;
vp = v/data.umean;

up = up + 0.01*randn(size(up)).*up;

% for t = 1:size(up,3)
%     up(:,:,t) = smoothn(up(:,:,t),'robust');
% end

h = fspecial('average',[5 5]);

for t = 1:size(up,3)
    up(:,:,t) = imfilter(up(:,:,t),h);
end

dt = 500;

Re = data.rho*data.umean*data.channelhgt/data.mu;

Po = data.rho*data.umean^2;

mask = u;
mask(mask~=0) = 1;
mask(:,1:3,:) = 0;
mask(:,end-2:end,:) = 0;
mask(1:3,:,:) = 0;
mask(end-2:end,:,:) = 0;

% se = strel('disk',3);
% mask = imerode(mask,se);
bcmask = permute(repmat(mask,[1 1 size(data.U,3)]),[2 1 3]);

solver = 'conservative';
%%
[p]=pressure_lines_9_2(xp',yp',permute(up,[2 1 3]),permute(vp,[2 1 3]),Re,dt,solver,bcmask(:,:,1));

%%
p = permute(p,[2 1 3])*Po;

pp_standard = p-repmat(p(size(data.x,1)/2,4,:),[size(data.x,1) size(data.x,2) 1]);