function crop_and_rotate(Data)
% This function rotates the image to obtain Cartesian flow and crops a
% rectangular region.
%
% The argument crop_window is a 1 x 4 vector containing [x0 y0 width height]
% where x0,y0 give the coordinates of the lower left corner of the region
% of interest.
%
% The argument angle is the angle of rotation in degrees.
%
% The two scalars frame_min and frame_max give the first and last frame in 
% the data sequence to process; if omitted the entire image sequence is 
% processed. Alternatively to load the entire image sequence
% to the end, frame_max may be set to Inf.
%
% Intensity scaling: none because cropping and rotation do not affect intensity

% extract relevant parameters from Data structure
read_directory  = Data.read_directory;
write_directory = Data.write_directory{Data.index};
image_format    = Data.image_format;
angle           = Data.crop_and_rotate.angle;
crop_window     = Data.crop_and_rotate.crop_window;
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;

direc_in = read_directory;
direc_1  = write_directory;

dlist = dir([direc_in,'*.',image_format]);        % list of images in the directory
if isinf(frame_max)
    frame_max = size(dlist,1);
end

% Rotation and cropping
for i=frame_min:frame_max
    fprintf('Image rotation and crop on Image %i\n',i);
    
    % read image
    I = imread(fullfile(direc_in,dlist(i).name));

    % image rotation (ir)
    % rotates image to align vessel axis with image axes, includes cropping
    I_ir = imrotate(I,angle,'bilinear'); 
    I_ir = imcrop(I_ir,crop_window);

    % write to file
    imwrite(I_ir,fullfile(direc_1,dlist(i).name));
end