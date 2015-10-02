function temporal_bandpass_old(Data)
% This function performs a high pass filter of the image set I along the 3rd
% dimension (time).  The function windows the data set with a Gaussian function to
% enforce periodicity and to minimize spectral leakage.
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
% Intensity scaling needs to be removed to be performed in a different
% filter
%
% Also needs to be fixed because a structure is used rather than a sequence
% of arguments

% extract relevant parameters from Data structure
read_directory  = Data.read_directory;
write_directory = Data.write_directory{Data.index};
image_format    = Data.image_format;
sigma_vect  = Data.temporal_bandpass.sigma_vect;
frame_min   = Data.frame_min;
frame_max   = Data.frame_max;

% Low- and high- pass frequency standard devation
lp_sigma=sigma_vect(1);
hp_sigma=sigma_vect(2);

% Maximum number of bytes that can be allotted to the data
% (leaving the *1e6 at the end converts the coefficient to MB)
max_mem=100*1e6;
% This is the list of images in the directory
image_list=dir([read_directory,'*.',image_format]);
% This returns an error if no images are found
if isempty(image_list)
    % This displays an error stating that no images were found
    error(['No images with the extension ',image_format,' were found in the specified directory.']);
end
% This checks whether the frame_min and frame_max arguments were supplied
% and if not sets the limits of the for loop to the entire sequence of
% images
if nargin==4
    % This is the minimum frame to load
    kk_min=1;
    % This is the maximum frame to load
    kk_max=length(image_list);
elseif nargin==6
    % This is the minimum frame to load
    kk_min=frame_min;
    % If the frame_max is set to infinity, then the last frame of the image
    % sequence is used, otherwise ii_max is set to frame_max
    if isinf(frame_max)
        % This sets the maximum frame to load to the length of the image
        % list
        kk_max=length(image_list);
    else
        % This sets the maximum frame to load to frame_max
        kk_max=frame_max;
    end
end
% This create the output directory if it does not exist
if not(exist(write_directory,'dir'))
    % This creates the output directory
    mkdir(write_directory);
end
% This loads the first image to extract the resolution and class type
I=imread([read_directory,image_list(kk_min).name]);
% This loads variable information about the image (including the image
% class)
image_info=whos('I');
% This reads the image class
image_class=image_info.class;
% This initializes the output image
I_Out=zeros(size(I));
% This converts the class of the output image to the same as the input
% images
eval(['I_Out=',image_class,'(I_Out);']);
% These are the image dimensions
[y_image_res,x_image_res]=size(I);
% This is the size of the windows to load (this assumes a square window)
win_res=round(sqrt(max_mem/(8*(kk_max-kk_min+1))));
% This is the number of windows to load in each dimension
x_win_num=ceil(x_image_res/win_res);
y_win_num=ceil(y_image_res/win_res);
% This is the total number of windows to load
win_num=x_win_num*y_win_num;
% This is the number of discrete samples in the time dimension
N=kk_max-kk_min+1;
% This is the standard deviation of the gaussian window
window_sigma=0.5;
% This is the windowing function to apply to the time series
W=exp(-0.5*(((0:N-1)-(N-1)/2)/(window_sigma*(N-1)/2)).^2);
% This permutes the windowing function to operate along the 3rd dimension
W=permute(W,[1,3,2]);

% This is a position vector to calculate a Gaussian mask (of the Fourier
% transformed image - this is separate from the Gaussian window function)
t_vector=2*(((1:N)-1)/(N-1)-0.5);
% This generates the Gaussian Fourier mask
H=1-1*exp(-(t_vector.^2)/(2*hp_sigma^2))-1*exp(-((abs(t_vector)-1).^2)/(2*lp_sigma^2));

% This permutes the Fourier filter to operate along the 3rd dimension
H=permute(H,[1,3,2]);
% This initializes the processed images minimum value
image_min=Inf;
% This initializes the processed images maximum value
image_max=-Inf;
% This iterates through the sections of the image to load
for n=1:win_num
    % These are the indices of the current window
    x_win_index=mod(n-1,x_win_num)+1;
    y_win_index=ceil(n/x_win_num);
    % This is the range of indices to read
    ii_min=(y_win_index-1)*win_res+1;
    ii_max=ii_min+win_res-1;
    jj_min=(x_win_index-1)*win_res+1;
    jj_max=jj_min+win_res-1;
    % This checks whether the maximum indices are greater than the size of
    % the image
    if ii_max>y_image_res
        ii_max=y_image_res;
    end
    if jj_max>x_image_res
        jj_max=x_image_res;
    end
    % This initializes the window matrix
    I_Win=zeros(ii_max-ii_min+1,jj_max-jj_min+1,kk_max-kk_min+1);
    % This iterates through the images
    for kk=kk_min:kk_max
        % This loads the current image
        I=imread([read_directory,image_list(kk).name]);
        % This extracts the roi
        I_Win(:,:,kk-kk_min+1)=double(I(ii_min:ii_max,jj_min:jj_max));
    end
    % This applies the windowing function to the current image sequence
    I_Win=bsxfun(@times,I_Win,W);
    % This converts the image to fourier space
    J=fftshift(fft(I_Win,[],3),3);
    % This performs the multiplication with the filter (convolution)
    K=bsxfun(@times,J,H);
    % This converts back to image space
    I_HP=ifft(ifftshift(K,3),[],3);
    % This removes residual imaginary parts
    I_HP=real(I_HP);
    % This rescales the image to invert the operation of the windowing;
    % this will likely cause images near the edge of the sequence to have
    % higher noise levels.
    I_HP=bsxfun(@rdivide,I_HP,W);
    
    % This is the minimum value of the current image sequence
    I_HP_Min=min(I_HP(:));
    % If this image sequence contains a minimum value less then the current
    % global value, the global value is reset
    if I_HP_Min<image_min
        % This resets the image minimum
        image_min=I_HP_Min;
    end
    % This is the maximum value of the current image sequence
    I_HP_Max=max(I_HP(:));
    % If this image sequence contains a maximum value greater then the current
    % global value, the global value is reset
    if I_HP_Max>image_max
        % This resets the image maximum
        image_max=I_HP_Max;
    end
    % This iterates through the list of output images saving them to the
    % output directory
    for m=1:N
        % If this is the first time through, then the output images are
        % initilized, if not then they are simply loaded from memory
        if n==1
            % This initializes the output matrix
            I_Out=zeros(y_image_res,x_image_res);
        else
            % This is the filename to read
            read_filename=[write_directory,image_list(m+kk_min-1).name];
            % This changes the output image extension to mat
            read_filename(end-2:end)='mat';
            % This loads in the current output matrix
            load(read_filename,'I_Out');
        end
        % This writes the current high pass filtered image window into the
        % saved output image
        I_Out(ii_min:ii_max,jj_min:jj_max)=I_HP(:,:,m);
        % This is the current filename to write out
        write_filename=[write_directory,image_list(m+kk_min-1).name];
        % This changes the output image extension to mat
        write_filename(end-2:end)='mat';
        % This writes the output image to memory
        save(write_filename,'I_Out');
    end
end
% This is the global intensity range
image_range=image_max-image_min;
% This iterates through the mat files rescaling them to fit between the
% global minimum and maximum values and saving the files as 16 bit tifs
for m=1:N
    % This is the filename to read
    read_filename=[write_directory,image_list(m+kk_min-1).name];
    % This changes the output image extension to mat
    read_filename(end-2:end)='mat';
    % This loads in the current output matrix
    load(read_filename,'I_Out');
    % This rescales I_Out by the global intensity range
    I_Out=(I_Out-image_min)/image_range;
    % This rescales I_Out to the 16 bit range
    I_Out=(2^16)*I_Out;
    % This converts I_Out to a uint16 image
    I_Out=uint16(I_Out);
    % This is the output filename
    write_filename=read_filename;
    % This changes the output image extension to tif
    write_filename(end-2:end)='tif';
    % This writes the output image to memory
    imwrite(I_Out,write_filename,'tif');
    % This deletes the current mat file
    delete(read_filename);
end