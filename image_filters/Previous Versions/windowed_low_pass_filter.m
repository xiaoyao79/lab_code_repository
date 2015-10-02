function windowed_low_pass_filter(read_directory,write_directory,image_format,lp_sigma,frame_min,frame_max);
% This function performs a low pass filter of the image I along the 3rd
% dimension.  The function windows the data set with a gaussian function to
% enforce periodicity and to minimize spectral leakage.  The argument
% lp_sigma is the standard deviation of the gaussian function centered at 
% zero frequency.  The frames are scaled between -1<=t<=1 so lp_sigma will
% be some fraction of 1 (ie lp_sigma=0.1).  The two scalars frame_min and 
% frame_max give the first and last frame in the data sequence to process; 
% these may be neglected to calculate the mean of the entire image sequence.
% Alternatively to load the entire image to the end, frame_max may be set
% to Inf.

% This is the maximum number of bytes that can be allotted to the data
% (leaving the *1e6 at the end converts the coefficient to MB)
max_mem=100*1e6;
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
if nargin==4;
    % This is the minimum frame to load
    kk_min=1;
    % This is the maximum frame to load
    kk_max=length(image_list);
elseif nargin==6;
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

% % This is the standard deviation of the gaussian low pass filter
% lp_sigma=0.1;

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