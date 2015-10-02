function image_sequence_filter_03;
% This function loads in a sequence of images in time and implements 1D
% filters on the pixels in time.  For standard data sets, the entire image
% sequence will not fit in the memory so only a window of each frame is
% loaded at one time.  These windows are then stacked together and the
% processing is completed on all pixels within the windows in parallel.
%
% Version 02 is similar to version 01 except that it is more functionalized
% to accomodate a wider array of operations including low-pass Fourier
% filtering, low-pass POD filtering, sequence median calculation, and
% sequence mean calculation.
%
% Version 03 uses nearly identical functions, but attempts to deal with
% scaling in a more appropriate manner, ie by saving as much information
% about the image sequence as possible given a type of compression.

% % This is the maximum number of bytes that can be allotted to the data
% max_mem=5e7;
% 
% % This is the directory to write the images to
% write_directory='/mnt/current_storage/Projects/Argonne/Argonne_Tests_072910/Calibration_02/Tubing_ReW0_11_5050Glycerin_Ga24_Test_01/Tomo_Reconstruction/Images/filtered_cam_1/';
% 
% write_directory='/home/rlafoy/Desktop/';
% 
% % This is the directory containing the images to read
% read_directory='/mnt/current_storage/Projects/Argonne/Argonne_Tests_072910/Calibration_02/Tubing_ReW0_11_5050Glycerin_Ga24_Test_01/Tomo_Reconstruction/Images/cam_1/';

% Good data sets from Argonne 2011 data:
%  /IRIS/Channel_Flow/cam_1/Calibration_01/PTFE_5UL_Test_07_Cam1/


read_directory='/home/rlafoy/Documents/AEThER/Argonne/2011_Data/cam_1/Calibration_01/PTFE_5UL_Test_07/';
write_directory='/home/rlafoy/Documents/AEThER/Argonne/2011_Data/cam_1/Calibration_01/PTFE_5UL_Test_07/filtered_images/';

tic;
windowed_band_pass_filter(read_directory,write_directory,'TIF',2,1001);

%I_Out=sequence_median(read_directory,'TIF',2,91);
disp(toc);
imagesc(I_Out);
axis equal;
axis image;
colormap('gray');



function I_Out=sequence_median(image_directory,image_format,frame_min,frame_max);
% This function calculates the median of a sequence of images in the
% directory given by the string image_directory (which must contain a slash
% at the end).  The string image_format is the extension of the images to
% search for in the directory.  The two scalars frame_min and frame_max
% give the first and last frame in the data sequence to process; these may
% be neglected to calculate the mean of the entire image sequence.
% Alternatively to load the entire image to the end, frame_max may be set
% to Inf.  The median is returned in the same format as the image format, ie 
% uint8 or uint16.

% This is the maximum number of bytes that can be allotted to the data
max_mem=1e9;
% This is the list of images in the directory
image_list=dir([image_directory,'*.',image_format]);
% This returns an error if no images are found
if isempty(image_list);
    % This displays an error stating that no images were found
    error(['No images with the extension ',image_format,' were found in the specified directory.']);
end;
% This checks whether the frame_min and frame_max arguments were supplied
% and if not sets the limits of the for loop to the entire sequence of
% images
if nargin==2;
    % This is the minimum frame to load
    kk_min=1;
    % This is the maximum frame to load
    kk_max=length(image_list);
elseif nargin==4;
    % This is the minimum frame to load
    kk_min=frame_min;
    % If the frame_max is set to infinity, then the last frame of the image
    % sequence is used, otherwise ii_max is set to frame_max
    if isinf(frame_max);
        % This sets the maximum frame to load to the length of the image
        % list
        kk_max=length(image_list);
    else;
        % This sets the maximum frame to load to frame_max
        kk_max=frame_max;
    end;
end;
% This loads the first image to extract the resolution and class type
I=imread([image_directory,image_list(kk_min).name]);
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
% These are the image resolutions
[y_image_res,x_image_res]=size(I);
% This is the size of the windows to load (this assumes a square window)
win_res=round(sqrt(max_mem/(8*(kk_max-kk_min+1))));
% This is the number of windows to load in each dimension
x_win_num=ceil(x_image_res/win_res);
y_win_num=ceil(y_image_res/win_res);
% This is the total number of windows to load
win_num=x_win_num*y_win_num;
% This iterates through the sections of the image to load
for n=1:win_num;
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
    if ii_max>y_image_res;
        ii_max=y_image_res;
    end;
    if jj_max>x_image_res;
        jj_max=x_image_res;
    end;
    % This initializes the window matrix
    I_Win=zeros(ii_max-ii_min+1,jj_max-jj_min+1,kk_max-kk_min+1);
    % This iterates through the images
    for kk=kk_min:kk_max;
        % This loads the current image
        I=imread([image_directory,image_list(kk).name]);
        % This extracts the roi
        I_Win(:,:,kk-kk_min+1)=I(ii_min:ii_max,jj_min:jj_max);
    end;
    % The next three commands utilize the MEX function fast_median wich
    % caclulates the median of a set in roughly O(n) versus O(n log n)
    % which is the performance speed of the matlab median command.  To use
    % the matlab median command, comment out the next three commands and
    % uncomment the fourth command.
    %
    % This permutes the I_Win matrix so that the fast_median function can
    % operate on the correct dimension of I_Win
    I_Win_Perm=permute(I_Win,[3,1,2]);
    % This calculates the median using the fast_median function
    I_Win_Median=fast_median(I_Win_Perm);
    % This reshapes the median image into the original dimensions
    I_Win_Median=reshape(I_Win_Median,size(I_Win,1),size(I_Win,2));
%     % This calculates the median of the time series of images
%     I_Win_Median=median(I_Win,3);
    % For values that are arithmetic means, this converts the output of the
    % median function to the same format as the input image
    eval(['I_Win_Median=',image_class,'(I_Win_Median);']);
    % This saves the median window to the output image
    I_Out(ii_min:ii_max,jj_min:jj_max)=I_Win_Median;
end;



function I_Out=sequence_mean(image_directory,image_format,frame_min,frame_max);
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
if isempty(image_list);
    % This displays an error stating that no images were found
    error(['No images with the extension ',image_format,' were found in the specified directory.']);
end;
% This checks whether the frame_min and frame_max arguments were supplied
% and if not sets the limits of the for loop to the entire sequence of
% images
if nargin==2;
    % This is the minimum frame to load
    ii_min=1;
    % This is the maximum frame to load
    ii_max=length(image_list);
elseif nargin==4;
    % This is the minimum frame to load
    ii_min=frame_min;
    % If the frame_max is set to infinity, then the last frame of the image
    % sequence is used, otherwise ii_max is set to frame_max
    if isinf(frame_max);
        % This sets the maximum frame to load to the length of the image
        % list
        ii_max=length(image_list);
    else;
        % This sets the maximum frame to load to frame_max
        ii_max=frame_max;
    end;
end;
% This iterates through the images to calculate the sequence mean
for ii=ii_min:ii_max;
    % This loads the first image in the sequence
    I=imread([image_directory,image_list(ii).name]);
    % If this is the first image to load in the sequence a summation matrix
    % is created and the image class is determined
    if ii==ii_min;
        % This initializes the summation matrix
        I_Sum=zeros(size(I));
        % This loads variable information about the image (including the image
        % class)
        image_info=whos('I');
        % This reads the image class
        image_class=image_info.class;
    end;
    % This adds the current image to the image sum
    I_Sum=I_Sum+double(I);
end;
% This divides the image sum by the total number of images loaded
I_Sum=I_Sum/(ii_max-ii_min+1);
% This converts the mean image to the format of the original images
eval(['I_Out=',image_class,'(I_Sum);']);



function windowed_low_pass_filter(read_directory,write_directory,image_format,frame_min,frame_max);
% This function performs a low pass filter of the image I along the 3rd
% dimension.  The function windows the data set with a gaussian function to
% enforce periodicity and to minimize spectral leakage.

% This is the maximum number of bytes that can be allotted to the data
max_mem=1e8;
% This is the list of images in the directory
image_list=dir([read_directory,'*.',image_format]);
% This returns an error if no images are found
if isempty(image_list);
    % This displays an error stating that no images were found
    error(['No images with the extension ',image_format,' were found in the specified directory.']);
end;
% This checks whether the frame_min and frame_max arguments were supplied
% and if not sets the limits of the for loop to the entire sequence of
% images
if nargin==3;
    % This is the minimum frame to load
    kk_min=1;
    % This is the maximum frame to load
    kk_max=length(image_list);
elseif nargin==5;
    % This is the minimum frame to load
    kk_min=frame_min;
    % If the frame_max is set to infinity, then the last frame of the image
    % sequence is used, otherwise ii_max is set to frame_max
    if isinf(frame_max);
        % This sets the maximum frame to load to the length of the image
        % list
        kk_max=length(image_list);
    else;
        % This sets the maximum frame to load to frame_max
        kk_max=frame_max;
    end;
end;
% This create the output directory if it does not exist
if not(exist(write_directory,'dir'));
    % This creates the output directory
    mkdir(write_directory);
end;
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
% These are the image resolutions
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
% This is the standard deviation of the gaussian low pass filter
lp_sigma=0.1;
% This is a position vector to calculate a gaussian mask (of the fourier
% transformed image - this is seperate than the gaussian window function)
t_vector=2*(((1:N)-1)/(N-1)-0.5);
% This generates the gaussian fourier mask
H=exp(-(t_vector.^2)/(2*lp_sigma^2));
% This permutes the fourier filter to operate along the 3rd dimension
H=permute(H,[1,3,2]);
% This iterates through the sections of the image to load
for n=1:win_num;
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
    if ii_max>y_image_res;
        ii_max=y_image_res;
    end;
    if jj_max>x_image_res;
        jj_max=x_image_res;
    end;
    % This initializes the window matrix
    I_Win=zeros(ii_max-ii_min+1,jj_max-jj_min+1,kk_max-kk_min+1);
    % This iterates through the images
    for kk=kk_min:kk_max;
        % This loads the current image
        I=imread([read_directory,image_list(kk).name]);
        % This extracts the roi
        I_Win(:,:,kk-kk_min+1)=double(I(ii_min:ii_max,jj_min:jj_max));
    end;
    % This applies the windowing function to the current image sequence
    I_Win=bsxfun(@times,I_Win,W);
    % This converts the image to fourier space
    J=fftshift(fft(I_Win,[],3),3);
    % This performs the multiplication with the filter (convolution)
    K=bsxfun(@times,J,H);
    % This converts back to image space
    I_LP=ifft(ifftshift(K,3),[],3);
    % This removes residual imaginary parts
    I_LP=real(I_LP);
    % This rescales the image to invert the operation of the windowing;
    % this will likely cause images near the edge of the sequence to have
    % higher noise levels.
    I_LP=bsxfun(@rdivide,I_LP,W);
    % This converts the image sequence back to the original format of the
    % imported images
    eval(['I_LP=',image_class,'(I_LP);']);
    % This iterates through the list of output images saving them to the
    % output directory
    for m=1:N;
        % If this is the first time through, then the output images are
        % initilized, if not then they are simply loaded from memory
        if n==1;
            % This initializes the output image
            I_Out=zeros(y_image_res,x_image_res);
            % This converts the image format to the same type as the input
            % image format
            eval(['I_Out=',image_class,'(I_Out);']);
        else;
            % This is the filename to read
            read_filename=[write_directory,image_list(m+kk_min-1).name];
            % This changes the output image extension to tif
            read_filename(end-2:end)='tif';
            % This reads in the current image from memory
            I_Out=imread(read_filename);
        end;
        % This writes the current low pass filtered image window into the
        % saved output image
        I_Out(ii_min:ii_max,jj_min:jj_max)=I_LP(:,:,m);
        % This is the current filename to write out
        write_filename=[write_directory,image_list(m+kk_min-1).name];
        % This changes the output image extension to tif
        write_filename(end-2:end)='tif';
        % This writes the output image to memory
        imwrite(I_Out,write_filename,'tif');
    end;
end;



function windowed_high_pass_filter(read_directory,write_directory,image_format,frame_min,frame_max);
% This function performs a high pass filter of the image I along the 3rd
% dimension.  The function windows the data set with a gaussian function to
% enforce periodicity and to minimize spectral leakage.

% This is the maximum number of bytes that can be allotted to the data
max_mem=1e8;
% This is the list of images in the directory
image_list=dir([read_directory,'*.',image_format]);
% This returns an error if no images are found
if isempty(image_list);
    % This displays an error stating that no images were found
    error(['No images with the extension ',image_format,' were found in the specified directory.']);
end;
% This checks whether the frame_min and frame_max arguments were supplied
% and if not sets the limits of the for loop to the entire sequence of
% images
if nargin==3;
    % This is the minimum frame to load
    kk_min=1;
    % This is the maximum frame to load
    kk_max=length(image_list);
elseif nargin==5;
    % This is the minimum frame to load
    kk_min=frame_min;
    % If the frame_max is set to infinity, then the last frame of the image
    % sequence is used, otherwise ii_max is set to frame_max
    if isinf(frame_max);
        % This sets the maximum frame to load to the length of the image
        % list
        kk_max=length(image_list);
    else;
        % This sets the maximum frame to load to frame_max
        kk_max=frame_max;
    end;
end;
% This create the output directory if it does not exist
if not(exist(write_directory,'dir'));
    % This creates the output directory
    mkdir(write_directory);
end;
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
% These are the image resolutions
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
% This is the standard deviation of the gaussian high pass filter
hp_sigma=0.1;
% This is a position vector to calculate a gaussian mask (of the fourier
% transformed image - this is seperate than the gaussian window function)
t_vector=2*(((1:N)-1)/(N-1)-0.5);
% This generates the gaussian fourier mask
H=1-exp(-(t_vector.^2)/(2*hp_sigma^2));
% This permutes the fourier filter to operate along the 3rd dimension
H=permute(H,[1,3,2]);

% This initializes the processed images minimum value
image_min=Inf;
% This initializes the processed images maximum value
image_max=-Inf;

% This iterates through the sections of the image to load
for n=1:win_num;
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
    if ii_max>y_image_res;
        ii_max=y_image_res;
    end;
    if jj_max>x_image_res;
        jj_max=x_image_res;
    end;
    % This initializes the window matrix
    I_Win=zeros(ii_max-ii_min+1,jj_max-jj_min+1,kk_max-kk_min+1);
    % This iterates through the images
    for kk=kk_min:kk_max;
        % This loads the current image
        I=imread([read_directory,image_list(kk).name]);
        % This extracts the roi
        I_Win(:,:,kk-kk_min+1)=double(I(ii_min:ii_max,jj_min:jj_max));
    end;
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
    if I_HP_Min<image_min;
        % This resets the image minimum
        image_min=I_HP_Min;
    end;
    % This is the maximum value of the current image sequence
    I_HP_Max=max(I_HP(:));
    % If this image sequence contains a maximum value greater then the current
    % global value, the global value is reset
    if I_HP_Max>image_max;
        % This resets the image maximum
        image_max=I_HP_Max;
    end;
    
    
    
%     % This converts the image sequence back to the original format of the
%     % imported images
%     eval(['I_HP=',image_class,'(I_HP);']);
    
    
    
    % This iterates through the list of output images saving them to the
    % output directory
    for m=1:N;
        % If this is the first time through, then the output images are
        % initilized, if not then they are simply loaded from memory
        if n==1;
            
%             % This initializes the output image
%             I_Out=zeros(y_image_res,x_image_res);
%             % This converts the image format to the same type as the input
%             % image format
%             eval(['I_Out=',image_class,'(I_Out);']);
            
            % This initializes the output matrix
            I_Out=zeros(y_image_res,x_image_res);
            
        else;
            
%             % This is the filename to read
%             read_filename=[write_directory,image_list(m+kk_min-1).name];
%             % This changes the output image extension to tif
%             read_filename(end-2:end)='tif';
%             % This reads in the current image from memory
%             I_Out=imread(read_filename);
            
            % This is the filename to read
            read_filename=[write_directory,image_list(m+kk_min-1).name];
            % This changes the output image extension to mat
            read_filename(end-2:end)='mat';
            % This loads in the current output matrix
            load(read_filename,'I_Out');
            
        end;
        
        
%         % This writes the current high pass filtered image window into the
%         % saved output image
%         I_Out(ii_min:ii_max,jj_min:jj_max)=I_HP(:,:,m);
%         % This is the current filename to write out
%         write_filename=[write_directory,image_list(m+kk_min-1).name];
%         % This changes the output image extension to tif
%         write_filename(end-2:end)='tif';
%         % This writes the output image to memory
%         imwrite(I_Out,write_filename,'tif');
        
        % This writes the current high pass filtered image window into the
        % saved output image
        I_Out(ii_min:ii_max,jj_min:jj_max)=I_HP(:,:,m);
        % This is the current filename to write out
        write_filename=[write_directory,image_list(m+kk_min-1).name];
        % This changes the output image extension to mat
        write_filename(end-2:end)='mat';
        % This writes the output image to memory
        save(write_filename,'I_Out');
        
    end;
end;
% This is the global intensity range
image_range=image_max-image_min;
% This iterates through the mat files rescaling them to fit between the
% global minimum and maximum values and saving the files as 16 bit tifs
for m=1:N;
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
end;




function windowed_band_pass_filter(read_directory,write_directory,image_format,frame_min,frame_max);
% This function performs a high pass filter of the image I along the 3rd
% dimension.  The function windows the data set with a gaussian function to
% enforce periodicity and to minimize spectral leakage.

% This is the maximum number of bytes that can be allotted to the data
max_mem=1e8;
% This is the list of images in the directory
image_list=dir([read_directory,'*.',image_format]);
% This returns an error if no images are found
if isempty(image_list);
    % This displays an error stating that no images were found
    error(['No images with the extension ',image_format,' were found in the specified directory.']);
end;
% This checks whether the frame_min and frame_max arguments were supplied
% and if not sets the limits of the for loop to the entire sequence of
% images
if nargin==3;
    % This is the minimum frame to load
    kk_min=1;
    % This is the maximum frame to load
    kk_max=length(image_list);
elseif nargin==5;
    % This is the minimum frame to load
    kk_min=frame_min;
    % If the frame_max is set to infinity, then the last frame of the image
    % sequence is used, otherwise ii_max is set to frame_max
    if isinf(frame_max);
        % This sets the maximum frame to load to the length of the image
        % list
        kk_max=length(image_list);
    else;
        % This sets the maximum frame to load to frame_max
        kk_max=frame_max;
    end;
end;
% This create the output directory if it does not exist
if not(exist(write_directory,'dir'));
    % This creates the output directory
    mkdir(write_directory);
end;
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
% These are the image resolutions
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


% This is the standard deviation of the gaussian high pass filter
hp_sigma=0.08;
% This is the standard deviation of the gaussian low pass filter
lp_sigma=0.02;
% This is a position vector to calculate a gaussian mask (of the fourier
% transformed image - this is seperate than the gaussian window function)
t_vector=2*(((1:N)-1)/(N-1)-0.5);
% This generates the gaussian fourier mask
H=1-1*exp(-(t_vector.^2)/(2*hp_sigma^2))-1*exp(-((abs(t_vector)-1).^2)/(2*lp_sigma^2));

% This permutes the fourier filter to operate along the 3rd dimension
H=permute(H,[1,3,2]);
% This initializes the processed images minimum value
image_min=Inf;
% This initializes the processed images maximum value
image_max=-Inf;
% This iterates through the sections of the image to load
for n=1:win_num;
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
    if ii_max>y_image_res;
        ii_max=y_image_res;
    end;
    if jj_max>x_image_res;
        jj_max=x_image_res;
    end;
    % This initializes the window matrix
    I_Win=zeros(ii_max-ii_min+1,jj_max-jj_min+1,kk_max-kk_min+1);
    % This iterates through the images
    for kk=kk_min:kk_max;
        % This loads the current image
        I=imread([read_directory,image_list(kk).name]);
        % This extracts the roi
        I_Win(:,:,kk-kk_min+1)=double(I(ii_min:ii_max,jj_min:jj_max));
    end;
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
    if I_HP_Min<image_min;
        % This resets the image minimum
        image_min=I_HP_Min;
    end;
    % This is the maximum value of the current image sequence
    I_HP_Max=max(I_HP(:));
    % If this image sequence contains a maximum value greater then the current
    % global value, the global value is reset
    if I_HP_Max>image_max;
        % This resets the image maximum
        image_max=I_HP_Max;
    end;
    % This iterates through the list of output images saving them to the
    % output directory
    for m=1:N;
        % If this is the first time through, then the output images are
        % initilized, if not then they are simply loaded from memory
        if n==1;
            % This initializes the output matrix
            I_Out=zeros(y_image_res,x_image_res);
        else;
            % This is the filename to read
            read_filename=[write_directory,image_list(m+kk_min-1).name];
            % This changes the output image extension to mat
            read_filename(end-2:end)='mat';
            % This loads in the current output matrix
            load(read_filename,'I_Out');
        end;
        % This writes the current high pass filtered image window into the
        % saved output image
        I_Out(ii_min:ii_max,jj_min:jj_max)=I_HP(:,:,m);
        % This is the current filename to write out
        write_filename=[write_directory,image_list(m+kk_min-1).name];
        % This changes the output image extension to mat
        write_filename(end-2:end)='mat';
        % This writes the output image to memory
        save(write_filename,'I_Out');
    end;
end;
% This is the global intensity range
image_range=image_max-image_min;
% This iterates through the mat files rescaling them to fit between the
% global minimum and maximum values and saving the files as 16 bit tifs
for m=1:N;
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
end;




function I_Out=low_pass_filter(I,sigma,image_class);
% This function performs a low pass filter of the image I along the 3rd
% dimension.

% This tiles the image to enforce periodicity
I=cat(3,I,flipdim(I,3));
% This converts the image to a double
I=double(I);
% This is the size of I
[yres,xres,zres]=size(I);
% This is a position vector to calculate a gaussian mask
x_vector=zeros(1,xres);
y_vector=zeros(1,yres);
z_vector=2*(((1:zres)-1)/(zres-1)-0.5);
% These are the position matricies
[X,Y,Z]=meshgrid(x_vector,y_vector,z_vector);
% This generates the gaussian mask
H=exp(-(Z.^2)/(2*sigma^2));
% This converts the image to fourier space
J=fftshift(fft(I,[],3),3);
% This performs the multiplication (convolution)
K=J.*H;
% This converts back to image space
I_Out=ifft(ifftshift(K,3),[],3);
% This removes residual imaginary parts
I_Out=real(I_Out);
% This extracts the first half of the image
I_Out=I_Out(:,:,1:size(I_Out,3)/2);
% This converts the image back to its orginal format
if strcmp(image_class,'uint8');
    I_Out=uint8(I_Out);
elseif strcmp(image_class,'uint16');
    I_Out=uint16(I_Out);
end;



function I_Out=pod_filter(I,mode_number);
% This function performs a singular value decomposition of I and then
% reconstructs the signal with a mode_number of modes

I=double(I);
x=zeros(size(I,3),size(I,1)*size(I,2));
for ii=1:size(I,3);
    I_Temp=I(:,:,ii);
    x(ii,:)=I_Temp(:)';
end;
[Ux,Sx,Vx]=svd(x);
x_modes=zeros(size(x));
for kk=1:mode_number;
    x_modes=x_modes+Ux(:,kk)*Sx(kk,kk)*(Vx(:,kk))';
end;

I_Out=zeros(size(I));
for ii=1:size(I,3);
    I_Out(:,:,ii)=reshape(x_modes(ii,:),size(I,1),size(I,2));
end;
I_Out=uint8(I_Out);

% I_Out=reshape(x_modes,size(I,1),size(I,2),size(I,3));
% I_Out=uint8(I_Out);

% I=double(I);
% I_Out=zeros(size(I));
% for ii=1:size(I,1);
%     for jj=1:size(I,2);
%         x=squeeze(I(ii,jj,:));
%         [Ux,Sx,Vx]=svd(x);
%         x_modes=zeros(size(x));
%         for kk=1:mode_number;
%             x_modes=x_modes+Ux(:,kk)*Sx(kk,kk)*(Vx(:,kk))';
%         end;
%         I_Out(ii,jj,:)=x_modes;
%     end;
% end;
% I_Out=uint8(I_Out);
