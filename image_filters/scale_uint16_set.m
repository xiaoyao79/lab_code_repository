function scale_uint16_set(Scaling)

% Function called by temporal_bandpass and spatial_bandpass to scale images
% in double format (saved as Matlab variables) to 16 bit integers with the
% image set global mean and maximum
%
% Also deletes the original .mat files.

image_max       = Scaling.image_max;
image_min       = Scaling.image_min;
N               = Scaling.N;
write_directory = Scaling.write_directory;
dlist           = Scaling.dlist;
frame_min       = Scaling.frame_min;
image_format    = Scaling.image_format;
image_name      = Scaling.image_name;

% Global intensity range
intensity_range=image_max-image_min;

% Loop through the mat files, rescaling between global minimum and maximum
% values and saving as 16 bit images in image_format
for m=1:N
    fprintf('Scaling image %i of %i\n',m,N)
    filename = [write_directory,dlist(m+frame_min-1).name];               % Generate file name based on original image list
    mat_filename = filename;
    mat_filename(end-size(image_format,2)+1:end) = 'mat';                   % Change the filename extension from image_format to 'mat'
    load(mat_filename,image_name);                                           % Open the variable
    
    eval(['I_Out = ',image_name,';']);
    
    % Scale the image to 16 bits with the global intensity range
    I_Out = (2^16)*(I_Out-image_min)/intensity_range;

    % Convert to uint16
    I_Out=uint16(I_Out);

    % Write the scaled 16 bit image file with the original filename
    imwrite(I_Out,filename);
    
    % Delete the mat file
    delete(mat_filename);
end