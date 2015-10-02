function submedian_inversion(Data)
% This function subtracts the image median from all pixels and inverts
% values below the median by computing the absolute value. This is intended
% to convert dark diffraction rings to bright positive signal
%%
% The two scalars frame_min and frame_max give the first and last frame in 
% the data sequence to process; if omitted the entire image sequence is 
% processed. Alternatively to load the entire image sequence
% to the end, frame_max may be set to Inf.
%
% Intensity scaling: Although the image is converted to double for
% inversion (uint8 and uint16 cannot be negative), no scaling is performed
% and the image is converted directly back to the original datatype after
% calculation of the absolute value

% extract relevant parameters from Data structure
read_directory  = Data.read_directory;
write_directory = Data.write_directory{Data.index};
image_format    = Data.image_format;
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;
    
dlist = dir([read_directory,'*.',image_format]);        % list of images in the directory
if isinf(frame_max)
    frame_max = size(dlist,1);
end

% Spatial highpass filter
for i=frame_min:frame_max
    fprintf('Boundary Inversion on Image %i\n',i);
    
    % Read image
    I = imread(fullfile(read_directory,dlist(i).name));
    
    image_class = class(I);

    I = double(I);
    
    % Perform median subtraction and inversion
    I2 = abs(I-median(I(:)));
    
    eval(['I2 =',image_class,'(I2);'])
    
    % write to file
    imwrite(I2,fullfile(write_directory,dlist(i).name));
end
