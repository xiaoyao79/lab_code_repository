function combineChannels(REDDIR, GREENDIR, BLUEDIR, IMBASERED, IMBASEGREEN, IMBASEBLUE, OUTDIR, OUTBASE, NDIGITS, EXT, START, STOP)
% combineChannels(REDDIR, GREENDIR, BLUEDIR, IMBASERED, IMBASEGREEN, IMBASEBLUE, OUTDIR, OUTBASE, NDIGITS, EXT, START, STOP)
% combines a series of separate color channel images (R G B) into a series of single color images.
%
% INPUTS
%   REDDIR = Path to directory containing red-channel images (string)
%
%   GREENDIR = Path to directory containing green-channel images (string) 
%
%   BLUEDIR = Path to directory containing blue-channel images (string) 
%
%   IMBASERED = Image base name for red-channel images (i.e., part of
%   the image name that preceeds image numbers) (string).
%
%   IMBASEGREEN = Image base name for green-channel images (string)
%
%   IMBASEBLUE = Image base name for blue-channel images (string)
%
%   OUTDIR = Path to directory in which color images will be saved (string). 
%
%   OUTBASE = Image base name for output color images (string)
%
%   NDIGITS = Number of digits in the image number suffix (integer)
%
%   EXT = Image extension, including 'dot' (example: EXT = '.tif') (string)
%
%   START = Number of the first image in the sequence (integer)
%
%   STOP = Number of the final image in the sequence (integer)
%
% OUTPUTS
%   None (images are saved to disk)
%
% SEE ALSO
%   filepaths

%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%

redPaths = filepaths(REDDIR, IMBASERED, NDIGITS, EXT, START, STOP); % Create paths to red images
greenPaths = filepaths(GREENDIR, IMBASEGREEN, NDIGITS, EXT, START, STOP); % Create paths to green images
bluePaths = filepaths(BLUEDIR, IMBASEBLUE, NDIGITS, EXT, START, STOP);  % Create paths to blue images
outPaths = filepaths(OUTDIR, OUTBASE, NDIGITS, EXT, START, STOP); % Create paths to color images

% Determine the number of images to be combined
nImages = size(redPaths, 1); 

% Determine image height and width (pixels)
[height width] = size(imread(redPaths(1, :))); 

% Combine the image channels
for k = 1:nImages
   
    colorImage = zeros(height, width, 3); % Initialize the color image
    redChannel = imread(redPaths(k, :)); % Read the red channel image
    greenChannel = imread(greenPaths(k, :)); % Read the green channel image
    blueChannel = imread(bluePaths(k, :)); % Read the blue channel image
    
    colorImage(:, :, 1) = redChannel; % Insert the red channel image into the first channel of the color image
    colorImage(:, :, 2) = greenChannel; % Insert the green channel image into the second channel of the color image
    colorImage(:, :, 3) = blueChannel; % Insert the blue channel image into the third channel of the color image
    
    imwrite(uint8(colorImage), outPaths(k, :)); % Write the color image to disk
    
end


end

%%%%%%%%%%%%%%%%%%%%
%%% END FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%


