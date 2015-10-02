function phase_retrieval(Data)
% This function performs phase retrieval as defined by Paganin et al 2004
%
% The argument tau is the tuning variable for the phase retrieval filter
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
frame_min       = Data.frame_min;
frame_max       = Data.frame_max;
tau             = Data.phase_retrieval.tau;
    
image_name = 'I2';

image_min = Inf;            % Initialize global minimum and maximum values of processed data set for final scaling
image_max = -Inf;

dlist = dir([read_directory,'*.',image_format]);        % list of images in the directory
if isinf(frame_max)
    frame_max = size(dlist,1);
end

% phase retrieval 
for i=frame_min:frame_max
    fprintf('Phase Retrieval of Image %i\n',i)
    
    % Read image
    I = double(imread(fullfile(read_directory,dlist(i).name)));
    
    % Perform phase retrieval (input and output are double format
    I2 = abs(pr(I,tau));
    
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

function I2 = pr(I,tau)
% Phase retrieval filter from Paganin 2004 (Eq 35) single image with thin object
% homogeneity assumption - does not use a regularization parameter as in
% Irvine 2008.  xi is the Fourier-space equivalent of the x- or y-coordinate
XX = size(I,2);
FFT_I = fft2(I);
frac = zeros(size(FFT_I));

for i=1:XX
    xi = i/XX;
    frac(:,i) = FFT_I(:,i)/(1+1j*xi*tau);
end

temp = ifft2(frac);
temp(temp<=0) = 1e-9;     % remove negative or zero values because the logarithm taken in the next step does not exist
I2 = real(-log(temp));

end

