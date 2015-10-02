function spatial_highpass(Data)
% This function performs a spatial high-pass filter using Gaussian windows
% on a set of images to remove illumination inhomogeneities.
%
% The argument sigma_vect is the standard deviation of the Gaussian 
% function, scaled to be between 0 and 1.
%
% The two scalars frame_min and frame_max give the first and last frame in 
% the data sequence to process; if omitted the entire image sequence is 
% processed. Alternatively to load the entire image sequence
% to the end, frame_max may be set to Inf.
%
% Intensity scaling: Images are written as individual mat files and scaled
% with scale_uint16_set

% extract relevant parameters from Data structure
read_directory  = Data.read_directory;
write_directory = Data.write_directory{Data.index};
image_format    = Data.image_format;
sigma           = Data.spatial_highpass.sigma;
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;
    
image_name = 'I2';

image_min = Inf;            % Initialize global minimum and maximum values of processed data set for final scaling
image_max = -Inf;

dlist = dir([read_directory,'*.',image_format]);        % list of images in the directory
if isinf(frame_max)
    frame_max = size(dlist,1);
end

% Spatial highpass filter
for i=frame_min:frame_max
    fprintf('Spatial Highpass Filtering Image %i\n',i);
    
    % Read image
    I = double(imread(fullfile(read_directory,dlist(i).name)));
    I = I(:,:,1);
    % Perform filter (input and output are double format)
    I2 = highpass_filter(I,sigma);
    
    % Update global image sequence maximum and minimum
    I_Min = min(I2(:));        % Minimum value of current block
    if image_min>I_Min           % If global minimum is greater than minimum of current block, change to minimum of current block 
        image_min = I_Min;
    end
    I_Max = max(I2(:));        % Maximum value of current block
    if image_max<I_Max           % If global maximum is greater than maximum of current block, change to maximum of current block 
        image_max = I_Max;
    end

    % Write variable to mat file
    filename = [write_directory,dlist(i).name];
    filename(end-size(image_format,2)+1:end) = 'mat';
    save(filename,image_name);
end

% Call post-processing scaling function to convert from doubles to uint16 images and delete the double files
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

%% Actual filter function
function I2 = highpass_filter(I,sigma)
% more_better_low_pass from Rod
% Tile the image to remove edge effects that arise during FFT
Itile=[fliplr(flipud(I)),flipud(I),fliplr(flipud(I));fliplr(I),I,fliplr(I);fliplr(flipud(I)),flipud(I),fliplr(flipud(I))]; %#ok<FLUDLR>

% Compute the image FFT
FFT_I = fftshift(fft2(Itile));

% Compute position vectors to create a grid for the Gaussian mask
[YY,XX]=size(Itile);
x_vector=2*(((1:XX)-1)/(XX-1)-0.5);
y_vector=2*(((1:YY)-1)/(YY-1)-0.5);
[X,Y]=meshgrid(x_vector,y_vector);

% Gaussian mask
H=exp(-(X.^2)/(2*sigma^2)-(Y.^2)/(2*sigma^2));

% Apply the Gaussian mask to the FFT of the image
FFT_I_masked=FFT_I.*H;

% Convert back to image space
Ilowpass=ifft2(ifftshift(FFT_I_masked));
I2=I-real(Ilowpass(YY/3+1:2*YY/3,XX/3+1:2*XX/3)); % extract real part of the center image

end