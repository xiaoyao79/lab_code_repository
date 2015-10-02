%%function highpass_butterworth(Data)

function discreate_wavelet_transform(Data)
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
levelnumber     = Data.levelnumber;
waveletname     = Data.waveletname;
alpha           = Data.alpha;
gbl_or_lvd      = Data.gbl_or_lvd;
s_or_h          = Data.s_or_h;
keep_approx     = Data.keep_approx;


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

%%
for i = frame_min:frame_max
    
    % Performs wavelet decomposition on image, calculates the approximate
    % coeffictient, and returns an image threshold for de-noising.
    [c,s] = wavedec2(ID(:,:,i),levelnumber,waveletname);
    A = appcoef2(c,s,waveletname,1);
    %thr = thselect(A,'rigrsure');
    [thr,nkeep] = wdcbm2(c,s,alpha,6*prod(s(1,:)));
    
    % Performs de-noising and reconstructs the now de-noised image from the
    % wavelet information.
    [xc,cxc,lxc,perf0,perf12] = wdencmp(gbl_or_lvd,c,s,waveletname,levelnumber,thr,s_or_h);
    X = uint8(xc);
    
    w  = fspecial('unsharp');
    Xs = double(X);
    Xf = uint8(imfilter(Xs,w,'replicate'));
    
    imwrite(Xf,fullfile(write_directory,dlist(i).name));
end