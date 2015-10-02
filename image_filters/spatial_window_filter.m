function spatial_window_filter(Data)
% This function performs temporal mean, median or minimum subtraction from an image
% set. It calls sequence_mean, sequence_median or sequence_minimum depending on the
% argument to_subtract and subtracts the resulting image from each image in
% the sequence.
%
% The two scalars frame_min and frame_max give the first and last frame in 
% the data sequence to process; if omitted the entire image sequence is 
% processed. Alternatively to load the entire image sequence
% to the end, frame_max may be set to Inf.
%
% Intensity scaling: none because the image subtracted is the same scaling
% as original images, giving a result with the same scaling as the original

% extract relevant parameters from Data structure
read_directory  = Data.read_directory;
write_directory = Data.write_directory{Data.index};
image_format    = Data.image_format;
to_subtract     = Data.background_subtraction.to_subtract;
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;

direc_in = read_directory;
direc_1  = write_directory;
dlist = dir([direc_in,'*.',image_format]);        % list of images in the directory
if isinf(frame_max)
    frame_max = size(dlist,1);
end

H = fspecial('sobel');

for i=frame_min:frame_max
    fprintf('Spatial Filter on Frame %i\n',i)
    I    = imread(fullfile(direc_in,dlist(i).name));
    Isub = imfilter(I,H);
    Isub = squeeze(Isub(:,:,1));
    imwrite(Isub,fullfile(direc_1,dlist(i).name));
end