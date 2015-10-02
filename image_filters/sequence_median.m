function I_Out=sequence_median(image_directory,image_format,frame_min,frame_max)
% This function calculates the median of a sequence of images in the
% directory given by the string image_directory (which must contain a slash
% at the end).  The string image_format is the extension of the images to
% search for in the directory.  The two scalars frame_min and frame_max
% give the first and last frame in the data sequence to process; these may
% be neglected to calculate the median of the entire image sequence.
% Alternatively to load the entire image to the end, frame_max may be set
% to Inf.  The median is returned in the same format as the image format, ie 
% uint8 or uint16.
%
% To use the fast_median function the mex files 'fast_median_mexa64' and 'fast_median'
% must be in the same directory as this function.

% This is the maximum number of bytes that can be allotted to the data
max_mem=1e9;
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
    kk_min=1;
    % This is the maximum frame to load
    kk_max=length(image_list);
elseif nargin==4
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
        I=imread([image_directory,image_list(kk).name]);
        % This extracts the roi
        I_Win(:,:,kk-kk_min+1)=I(ii_min:ii_max,jj_min:jj_max);
    end
    % The next three commands utilize the MEX function fast_median wich
    % caclulates the median of a set in roughly O(n) versus O(n log n)
    % which is the performance speed of the matlab median command.  To use
    % the matlab median command, comment out the next three commands and
    % uncomment the fourth command.
    %
%     % This permutes the I_Win matrix so that the fast_median function can
%     % operate on the correct dimension of I_Win
%     I_Win_Perm=permute(I_Win,[3,1,2]);
%     % This calculates the median using the fast_median function
%     I_Win_Median=fast_median(I_Win_Perm);
%     % This reshapes the median image into the original dimensions
%     I_Win_Median=reshape(I_Win_Median,size(I_Win,1),size(I_Win,2));
    % This calculates the median of the time series of images
    I_Win_Median=median(I_Win,3);
    % For values that are arithmetic means, this converts the output of the
    % median function to the same format as the input image
    eval(['I_Win_Median=',image_class,'(I_Win_Median);']);
    % This saves the median window to the output image
    I_Out(ii_min:ii_max,jj_min:jj_max)=I_Win_Median;
end