%%function highpass_butterworth(Data)

function median_window_sub(Data)
% This function performs Butterworth filtering on an image set I along the 3rd
% dimension (time) through the frequency domain (Hz).  The function transforms
% the intensity information for each frame into the frequency domain.  A
% 3rd order Butterworth filter at a normalized cut-off frequency of
% is then applied to each pixel along the frequency domain.  The filtered
% image set is then transformed back into the time domain, the real
% elements of the double are taken and the image is reintensified and made
% in to a uint8 data set again.
%
% The argument n_order is the order approximation that the filter work
% toward.  For example, if the filter is set to a 3rd order approximation
% the filter will run through the approximation algorithm until it reaches
% 3rd order.  The second argument, cutoff_freq, is the normalized cuttoff
% for the data set.  This value must be within 0 and 1, with a typical
% start value of 0.33.  The final argument ftype give the type of filter
% you wish to set.  These can either be 'low','high' or 'stop'.
%
% The two scalars frame_min and frame_max give the first and last frame in
% the data sequence to process; if omitted the entire image sequence is
% processed. Alternatively to load the entire image sequence
% to the end, frame_max may be set to Inf.
%
% Intensity scaling has been moved to a subfunction (scale_uint16_set
% usable by other filters such as spatial_highpass
%
% File created for use with filter_main on 12/8/2011 by Brett Meyers


read_directory  = Data.read_directory;
write_directory = Data.write_directory{Data.index};
image_format    = Data.image_format;
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;
windowsize      = Data.windowsize;
to_subtract     = Data.median_window_sub.to_subtract;


image_name = 'I2';

image_min = Inf;            % Initialize global minimum and maximum values of processed data set for final scaling
image_max = -Inf;

dlist = dir([read_directory,'*.',image_format]);        % list of images in the directory
if isinf(frame_max)
    frame_max = size(dlist,1);
end

fprintf('Median Window Subtraction Filtering In Progress');

for t=frame_min-windowsize:frame_max
    fprintf('Reading frame %i\n',t);
    
    % Read image
    I(:,:,:,t)  = imread(fullfile(read_directory,dlist(t).name));
    ID(:,:,t)   = squeeze(I(:,:,1,t));
    %     ID(:,:,t) = double(I(:,:,t));
end
%% Creates an identity matrix
time_mat=diag(ones(frame_max,1),0);

% for-loop adds "ones" to the right and left of each diagonal of the identity matrix
for k = 1:windowsize
    time_mat=time_mat+diag(ones(frame_max-(k),1), k)+...
        diag(ones(frame_max-(k),1),-k);
end

% Creates a matrix of values which correspond to each image in the image set
value_mat = cumsum(ones(frame_max),2);

% Sets any element in value_mat equal to zero if element in time_mat is zero
value_mat(time_mat==0)=0;

% Generates values that correspond to value_mat size dimensions
[row,col]=size(value_mat);

% Generates values that correspond to the original image dimensions
[dim1,dim2,dim3]=size(ID);

% Generates initial value matrices for mean, median and minimum filter arrays
% filtmean = zeros([dim1,dim2,row],'uint8');
filtmed  = zeros([dim1,dim2,row],'uint8');
% filtmin  = zeros([dim1,dim2,row],'uint8');

%% for-loop creates the filter image arrays and saves them to .tif files
for n=frame_min:(frame_max-windowsize)
    % Creates an array of numbers to chose which original image array values are used for filtering
    s=nonzeros(value_mat(n,:));
    s=s';
    % Image filtering step
    switch to_subtract
        case 'mean'
            filtframe   = uint8(mean(ID(:,:,s),3));
            filtmean    = ID(:,:,value_mat(n,n)) - uint8(mean(ID(:,:,s),3));
            scaleFrame  = max(max(filtmean));
            filtmean    = double(filtmean);
            scaleFrame  = double(scaleFrame);
            filtmean    = 255*(filtmean./scaleFrame);
            filtmean    = uint8(filtmean);
            imwrite(filtmean,fullfile(write_directory,dlist(n).name));
            imwrite(filtframe,fullfile([write_directory,'filter ',dlist(n).name]));

        case 'median'
            filtmed    = ID(:,:,value_mat(n,n)) - uint8(median(ID(:,:,s),3));
            scaleFrame = max(max(filtmed));
            filtmed    = double(filtmed);
            scaleFrame = double(scaleFrame);
            filtmed    = 255*(filtmed./scaleFrame);
            filtmed    = uint8(filtmed);
            imwrite(filtmed,fullfile(write_directory,dlist(n).name));
        case 'minimum'
            filtmin    = ID(:,:,value_mat(n,n)) - uint8(min(ID(:,:,s),[],3));
            scaleFrame = max(max(filtmin));
            filtmin    = double(filtmin);
            scaleFrame = double(scaleFrame);
            filtmin    = 255*(filtmin./scaleFrame);
            filtmin    = uint8(filtmin);
            imwrite(filtmin,fullfile(write_directory,dlist(n).name));
    end
    
end

end