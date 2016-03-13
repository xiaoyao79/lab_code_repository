function [data] = cart_pol_grid_gen(data,info)
% Subfunction for defining the cartesian and polar grids to be used in
% image analysis including (1) image unwrapping and (2) field
% reconstruction

% Select corners of the Color Doppler window to determine the polar grid
% zero-location
temp = figure(1); 
imagesc(uint8(data.imcdopp(:,:,:,1)));
set(gcf,'position',[100 100 0.6*info.s2 0.6*info.s1]);
fprintf('Select region edges (4 points needed) to unwrap window \r')
[xloc,yloc] = getpts(temp);

% Find the location of the polar grid zero-location
% Identify left-most points
xsort = sort(xloc,'ascend');
for i = 1:2
    r = find(xloc == xsort(i));
    xleft(i) = xloc(r); 
    yleft(i) = yloc(r);
end
p1 = polyfit(xleft,yleft,1);

% Identify right-most points
xsort = sort(xloc,'descend');
for i = 1:2
    r = find(xloc == xsort(i));
    xright(i) = xloc(r); 
    yright(i) = yloc(r);
end
p2 = polyfit(xright,yright,1);

% Find intersections of the lines. This will identify the start of the
% sweep angle region for extraction of the color doppler signal subregion
% in the image sequence
x_intersect = round(fzero(@(x) polyval(p1-p2,x),3));
y_intersect = round(polyval(p1,x_intersect));

% Generate the Cartesian Grid of the image (relative to sweep window
% center)
[data.X,data.Y] = meshgrid(linspace(1-x_intersect,info.s2-x_intersect,info.s2),...
    linspace(1-y_intersect,info.s1,info.s1));

% Generate Polar Grid
[data.THETA,data.RHO] = cart2pol(data.X,data.Y);
data.Theta = radtodeg(data.THETA);
data.Rho   = data.RHO;

% Store point locations (used later in extracting regions)
data.xloc  = [xleft';xright'];
data.yloc  = [yleft';yright'];

% try
%     Theta = Theta(end-data.s1+1:end,:);
%     Rho   = RHO(end-data.s1+1:end,:);
%     
%     THETA = THETA(end-data.s1+1:end,:);
%     RHO   = RHO(end-data.s1+1:end,:);
% catch
%     temp = zeros(size(data.imcdopp(:,:,1)));
%     temp(end-size(Theta,1)+1:end,:) = Theta;
%     
%     Theta = temp;
%     
%     temp = zeros(size(data.imcdopp(:,:,1)));
%     temp(end-size(RHO,1)+1:end,:) = RHO;
%     
%     Rho   = temp;
%     
%     temp = zeros(size(data.imcdopp(:,:,1)));
%     temp(end-size(THETA,1)+1:end,:) = THETA;
%     
%     THETA = temp;
%     
%     RHO   = Rho;
% end
% 
% data.Theta = Theta;
% data.Rho   = Rho;
% data.THETA = THETA;
% data.RHO   = RHO;
% data.xloc  = xloc;
% data.yloc  = yloc;