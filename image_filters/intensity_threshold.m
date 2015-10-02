function intensity_threshold(Data)
% This function thresholds the data intensity range to enhance contrast of
% middle values. Only supports images with class uint8 or uint16.
%
% The argument threshold is a 1 x 2 vector containing [lower_threshold
% upper_threshold] where each threshold is a fraction of the total dynamic
% range.
%
% The two scalars frame_min and frame_max give the first and last frame in 
% the data sequence to process; if omitted the entire image sequence is 
% processed. Alternatively to load the entire image sequence
% to the end, frame_max may be set to Inf.
%
% Intensity scaling: maximum and minimum thresholding

% extract relevant parameters from Data structure
read_directory  = Data.read_directory;
write_directory = Data.write_directory{Data.index};
image_format    = Data.image_format;
lower_threshold = Data.intensity_threshold.threshold(1);
upper_threshold = Data.intensity_threshold.threshold(2);
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;

dlist = dir([read_directory,'*.',image_format]);        % list of images in the directory
if isinf(frame_max)
    frame_max = size(dlist,1);
end

% Read sample image to get image class
Itemp = imread(fullfile(read_directory,dlist(1).name));
image_class = class(Itemp);
image;

if strcmp(image_class,'uint8')
    range = 2^8;
else if strcmp(image_class,'uint16')
        range = 2^16;
    else error('Image class %s not supported\n',image_class)
    end
end

% Set new maximum and minimum intensity based on thresholds
maxI = upper_threshold*(range-1);
minI = lower_threshold*(range-1);

% Thresholding
for i=frame_min:frame_max
    fprintf('Thresholding Image %i\n',i);
    
    % Read image
    I = double(imread(fullfile(read_directory,dlist(i).name)));
    
    % Image thresholding
    I2 = (range-1)*(I-minI)/(maxI-minI);
    
    % Ensure that values outside thresholds are written as new maximum and
    % minimum
    I2(I2<0) = 0;
    I2(I2>range-1) = range-1;
    
    % Conversion to original data type
    eval(['I2 =',image_class,'(I2);'])

    % write to file
    imwrite(I2,fullfile(write_directory,dlist(i).name));
end

end