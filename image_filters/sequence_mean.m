function I_Out=sequence_mean(image_directory,image_format,frame_min,frame_max) %#ok<STOUT> This is set in the final eval statement
% This function calculates the mean of a sequence of images in the
% directory given by the string image_directory (which must contain a slash
% at the end).  The string image_format is the extension of the images to
% search for in the directory.  The two scalars frame_min and frame_max
% give the first and last frame in the data sequence to process; these may
% be neglected to calculate the mean of the entire image sequence.
% Alternatively to load the entire image to the end, frame_max may be set
% to Inf.  The mean is returned in the same format as the image format, ie 
% uint8 or uint16.

% This is the list of images in the directory
image_list=dir([image_directory,'*.',image_format]);
% This returns an error if no images are found
if isempty(image_list)
    % This displays an error stating that no images were found
    error(['No images with the extension ',image_format,' were found in the specified directory.']);
end
% This checks whether the frame_min and frame_max arguments were supplied
% and if not sets the limits of the for loop to the entire sequence of
% images
if nargin==2
    % This is the minimum frame to load
    ii_min=1;
    % This is the maximum frame to load
    ii_max=length(image_list);
elseif nargin==4
    % This is the minimum frame to load
    ii_min=frame_min;
    % If the frame_max is set to infinity, then the last frame of the image
    % sequence is used, otherwise ii_max is set to frame_max
    if isinf(frame_max)
        % This sets the maximum frame to load to the length of the image
        % list
        ii_max=length(image_list);
    else
        % This sets the maximum frame to load to frame_max
        ii_max=frame_max;
    end
end
% This iterates through the images to calculate the sequence mean
for ii=ii_min:ii_max
    % This loads the first image in the sequence
    I=imread([image_directory,image_list(ii).name]);
    % If this is the first image to load in the sequence a summation matrix
    % is created and the image class is determined
    if ii==ii_min
        % This initializes the summation matrix
        I_Sum=zeros(size(I));
        % This loads variable information about the image (including the image
        % class)
        image_info=whos('I');
        % This reads the image class
        image_class=image_info.class;
    end
    % This adds the current image to the image sum
    I_Sum=I_Sum+double(I);
end
% This divides the image sum by the total number of images loaded
I_Sum=I_Sum/(ii_max-ii_min+1); %#ok<NASGU> This is used in the final eval statement
% This converts the mean image to the format of the original images
eval(['I_Out=',image_class,'(I_Sum);']);