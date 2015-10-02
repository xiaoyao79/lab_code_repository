function temporal_bandpass(Data)
% This function performs a high pass filter of the image set I along the 3rd
% dimension (time).  The function windows the data set with a Gaussian function to
% enforce periodicity and to minimize spectral leakage. Because this is a
% temporal function operating on single pixels throughout the entire data
% set, block processing is implemented to avoid the need to load the entire
% image set into memory (which is likely to overload system memory).
%
% The argument sigma_vect is a 1 x 2 vector containing [lp_sigma,hp_sigma]
% or the standard deviation of the gaussian function centered at  zero 
% frequency and at the frame_rate frequency.  The frames are scaled between
% -1<=t<=1 so lp_sigma and hp_sigma will be some fraction of 1 (ie 
% lp_sigma=0.02 and hp_sigma=0.08). Set either value to 0 to obtain purely
% a high-pass or low-pass filter.
%
% The two scalars frame_min and frame_max give the first and last frame in 
% the data sequence to process; if omitted the entire image sequence is 
% processed. Alternatively to load the entire image sequence
% to the end, frame_max may be set to Inf.
%
% This function was previously called windowed_band_pass_filter, renamed to
% differentiate between this (temporal) and the spatial high pass filter
% and modified to accept the structure argument.
%
% Intensity scaling has been moved to a subfunction (scale_uint16_set 
% usable by other filters such as spatial_highpass
%
% Eventually the value of sigma for the time signal Gaussian window may
% become an argument - currently hard-coded to be 0.5

% extract relevant parameters from Data structure
read_directory  = Data.read_directory;
write_directory = Data.write_directory{Data.index};
image_format    = Data.image_format;
lp_sigma        = Data.temporal_bandpass.sigma_vect(1);
hp_sigma        = Data.temporal_bandpass.sigma_vect(2);
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;


%% Minor calculation and initialization
max_mem = 100*1e6;          % Memory allocation for filter - default 100MB

image_min = Inf;            % Initialize global minimum and maximum values of processed data set for final scaling
image_max = -Inf;

image_name = 'I_Out';

%% Directory input
dlist=dir([read_directory,'*.',image_format]);
if isempty(dlist)
    error(['No images with the extension ',image_format,' were found in the specified directory.']);
end
if isinf(frame_max)
    frame_max = size(dlist,1);
end
N = frame_max-frame_min+1;  % Total number of images

%% Load single image to extract dimensions
Itemp = imread(fullfile(read_directory,dlist(frame_min).name));
[height,width] = size(Itemp);

%% Parameters for block computation to avoid memory overload
block_size  = round(sqrt(max_mem/(8*(frame_max-frame_min+1))));     % Size of blocks to load (assuming square window) - determined by size of image set and memory allocation
x_block_num = ceil(width/block_size);                               % Number of blocks to load in each dimension
y_block_num = ceil(height/block_size);
block_num   = x_block_num*y_block_num;                              % Total number of blocks to load

%% Filter parameters
% Gaussian window on time series signal to enforce periodicity of the data
window_sigma = 0.5;                                                         % Standard deviation
W = exp(-0.5*(((0:N-1)-(N-1)/2)/(window_sigma*(N-1)/2)).^2);                % Windowing function
W = permute(W,[1,3,2]);                                                     % Permutation so that the function operates in the 3rd dimension (time)

% Gaussian mask on Fourier transformed image
t_vector = 2*(((1:N)-1)/(N-1)-0.5);                                         % Position vector for mask
H = 1-1*exp(-(t_vector.^2)/(2*hp_sigma^2))-1*exp(-((abs(t_vector)-1).^2)/(2*lp_sigma^2));   % Gaussian mask computed from position vector and high- and low-pass standard deviations
H = permute(H,[1,3,2]);                                                     % Permutation so that the function operates in the 3rd dimension (time)

%% Filter application
% Iterate through the image blocks (defined in parameters for block
% computation)
for n=1:block_num
    fprintf('Filtering block %i of %i\n',n,block_num);
    
    % Overall indices of the current block
    x_block_index = mod(n-1,x_block_num)+1;
    y_block_index = ceil(n/x_block_num);
    
    % Range of indices to read (pixel locations within the current block)
    ii_min = (y_block_index-1)*block_size+1;
    ii_max = ii_min+block_size-1;
    jj_min = (x_block_index-1)*block_size+1;
    jj_max = jj_min+block_size-1;
    
    % Limit maximum indices to the image dimensions
    if ii_max>height
        ii_max = height;
    end
    if jj_max>width
        jj_max = width;
    end
    
    % Initialize the block
    I_Win = zeros(ii_max-ii_min+1,jj_max-jj_min+1,N);
    
    % Load the block for all images
    for i=frame_min:frame_max
        I = imread([read_directory,dlist(i).name]);                         % Load entire image
        I_Win(:,:,i-frame_min+1) = double(I(ii_min:ii_max,jj_min:jj_max));  % Extract block
    end
    
    %% Filter
    I_Win = bsxfun(@times,I_Win,W);                         % Temporal windowing of the current block = why bsxfun?
    J = fftshift(fft(I_Win,[],3),3);                        % Transform to Fourier space
    K = bsxfun(@times,J,H);                                 % Convolve signal with Fourier mask
    I_HP = ifft(ifftshift(K,3),[],3);                       % Inverse transform back to image space
    I_HP = real(I_HP);                                      % Remove residual imaginary parts
    
    % Rescale block to invert the time series windowing operation.
    % This will likely cause images near the edge of the sequence to have
    % higher noise levels.
    I_HP = bsxfun(@rdivide,I_HP,W);       
    
    %% Update global image sequence maximum and minimum
    I_HP_Min = min(I_HP(:));        % Minimum value of current block
    if image_min>I_HP_Min           % If global minimum is greater than minimum of current block, change to minimum of current block 
        image_min = I_HP_Min;
    end
    I_HP_Max = max(I_HP(:));        % Maximum value of current block
    if image_max<I_HP_Max           % If global maximum is greater than maximum of current block, change to maximum of current block 
        image_max = I_HP_Max;
    end
    
    %% Update (write) output images as variable in mat files to retain double class
    % Iterate through the list of images
    for m=1:N
        % Get .mat file name containing image variable
        filename = [write_directory,dlist(m+frame_min-1).name];           % Generate file name based on original image list
        filename(end-size(image_format,2)+1:end) = 'mat';                 % Change the filename extension from image_format to 'mat'
        
        if n==1                      
            I_Out = zeros(size(Itemp));                                   % Initialize the output variable if this is the first loop
        else                       
            load(filename,image_name);                                    % Otherwise open the existing variable                                     
        end

        I_Out(ii_min:ii_max,jj_min:jj_max) = I_HP(:,:,m);                 % Add the filtered block to the output variable

        save(filename,image_name);                                        % Save the current block to file
    end
end

% Call post-processing scaling function to convert from doubles to uint16 images and delete the double files
Scaling.image_max = image_max;
Scaling.image_min = image_min;
Scaling.N = N;
Scaling.write_directory = write_directory;
Scaling.dlist = dlist;
Scaling.frame_min = frame_min;
Scaling.image_format = image_format;
Scaling.image_name = image_name;

scale_uint16_set(Scaling)

end