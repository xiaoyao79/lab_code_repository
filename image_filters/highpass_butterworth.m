%%function highpass_butterworth(Data)

function highpass_butterworth(Data)
% This function performs a highpass filter of the image set I along the 3rd
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
% 3rd order.  The other argument, cutoff_freq, is the normalized cuttoff 
% for the data set.  This value must be within 0 and 1, with a typical
% start value of 0.33.
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
n_order         = Data.n_order;
cutoff_freq     = Data.cutoff_freq;
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;

image_name = 'I2';

image_min = Inf;            % Initialize global minimum and maximum values of processed data set for final scaling
image_max = -Inf;

dlist = dir([read_directory,'*.',image_format]);        % list of images in the directory
if isinf(frame_max)
    frame_max = size(dlist,1);
end

fprintf('Highpass Butterworth Filtering In Progress');

for t=frame_min:frame_max
    fprintf('Reading frame %i\n',t);
    
    % Read image
    I(:,:,t)  = imread(fullfile(read_directory,dlist(t).name));
    ID(:,:,t) = double(I(:,:,t));
end

%%  Fast Fourier Transfor of Image set
for t = 1:size(ID,3)
    FFT(:,:,t)=fftshift(fft2(ID(:,:,t),2*size(ID,1),2*size(ID,2)));
end

%%  Butterworth Filter Design
[b,a]=butter(3,1/3,'high');

%%  Filtering Stage
for j = 1:size(FFT,1)
    for i = 1:size(FFT,2)
        Filt_FFT(j,i,:)=filtfilt(b,a,FFT(j,i,:));
    end
end

%%  Inverse Fourier Transfrom
for t = 1:size(FFT,3)
    fprintf('Returning Frame %i\n',t)
    IFFT(:,:,t)=ifft2(ifftshift(Filt_FFT(:,:,t)));
    I2(:,:,t)=uint8(real(IFFT(1:size(normal,1),1:size(normal,2),t)));
    
%%  Write variable to mat file
    filename = [write_directory,dlist(t).name];
    filename(end-size(image_format,2)+1:end) = 'mat';
    save(filename,image_name);
end

%%  Filtering Preview
for t = 1:size(I2,3)
    figure(1)
    title('%i/n',t);
    imagesc(filt_renorm(:,:,t))
    map=colormap(gray)
end

%%  Update global image sequence maximum and minimum
I_Min = min(I2(:));        % Minimum value of current block
if image_min>I_Min           % If global minimum is greater than minimum of current block, change to minimum of current block
    image_min = I_Min;
end
I_Max = max(I2(:));        % Maximum value of current block
if image_max<I_Max           % If global maximum is greater than maximum of current block, change to maximum of current block
    image_max = I_Max;
end

Scaling.image_max       = image_max;
Scaling.image_min       = image_min;
Scaling.N               = frame_max-frame_min+1;
Scaling.write_directory = write_directory;
Scaling.dlist           = dlist;
Scaling.frame_min       = frame_min;
Scaling.image_format    = image_format;
Scaling.image_name      = image_name;

scale_uint16_set(Scaling)
end