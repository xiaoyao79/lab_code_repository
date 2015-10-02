function hessianRun(IMDIR, IMBASE, OUTDIR, OUTBASE, NDIGITS, EXT, START, STOP, GAUSSIANSCALES, KERNELRADIUS, RESCALE, SCALINGCOEFFICIENTS, NPROCESSORS)
% hessianRun(IMDIR, IMBASE, OUTDIR, OUTBASE, NDIGITS, EXT, START, STOP, GAUSSIANSCALES, KERNELRADIUS, RESCALE, SCALINGCOEFFICIENTS, NPROCESSORS)
% Calculates the 3-D hessian eigen-images of a series of images.
%
% KNOWN ISSUES
%   Parallel processing in this function is incompatable with arbitrary
%   derivative convolution kernel lengths right now. Therefore, the
%   convolution kernel radius is forced to 4 when NPROCESSORS > 1.
% 
% INPUTS
%   IMDIR = Path to directory in which images are located (string).
%
%   IMBASE = Image base name (i.e., part of the image name that preceeds
%   image numbers) (string).
%
%   OUTDIR = Directory underneath which the output directories will be
%   created. Processed images are saved to these subdirectories. (string). 
%
%   OUTBASE = Base name of output images (i.e., part of the image name that preceeds
%   image numbers) (string).
%
%   NDIGITS = Number of digits in the image number suffix (integer)
%
%   EXT = Image extension, including 'dot' (example: EXT = '.tif') (string)
%
%   START = Number of the first image in the sequence (integer)
%
%   STOP = Number of the final image in the sequence (integer)
%
%   GAUSSIANSCALES = [1 x n] vector specifying the standard deviation
%   (in pixels) of the gaussian kernels used to calculate image
%   derivatives, where n = number of dimensions in the image (n = 2 or 3).
%   If GAUSSIANSCALES is not specified, it defaults to 1 pixel in every
%   direction, i.e., GAUSSIANSCALES = [ 1 1 1 ]. 
%
%  KERNELRADIUS = Radius (in pixels) of convolution kernel (integer). The
%  kernel radius should be an an even integer so that the kernel is symmetric
%  about its anchor point; therefore, if an odd or non-integer kernel
%  radius is specified, the kernel radius is forced to the next greatest
%  even integer. 
% 
%  RESCALE = Binary flag that specifies whether or not to scale the
%   Z-derivatives by a constant value such that the order of the magnitude
%   of the  Z derivatives is similar to those of the X- and Y- derivatives. 
%
%  SCALINGCOEFFICIENTS = Vector containing values of forced ratios of
%   the mean magnitude of the Z- intensity derivative to those of the
%   X and Y intensity derivatives. 
%
%   NPROCESSORS = Number of processors to use in parallel (integer).
%   
% OUTPUTS
%   None (images are saved to disk)
%
% SEE ALSO
%   hessian3, checkHessianJobMemory, filepaths, combineChannels

%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% % SET DEFAULTS %%
%%%%%%%%%%%%

% Default to 1 processor
if nargin < 13
    NPROCESSORS = 1;
end

% Default to scaling coefficient of 1;
if nargin < 12
    SCALINGCOEFFICIENTS = 1;
end

% Default to no rescaling
if nargin < 11
    RESCALE = 0;
end

% Default to a kernel radius of 4 pixels
if nargin < 10
    KERNELRADIUS = 4;
end

% Default to gaussian kernel scales of 1 pixel in every direction
if nargin < 9
    GAUSSIANSCALES = [1 1 1];
end

% Create paths to raw images
paths = filepaths(IMDIR, IMBASE, NDIGITS, EXT, START, STOP); 

% Calculate the number of images
nImages = STOP - START + 1; 

% Determine size of images
[height width nChannels] = size(imread(paths(1, :))); 

% Create output basenames for each color
outbaseRed = [OUTBASE 'red_'];
outbaseGreen = [OUTBASE 'green_'];
outbaseBlue = [OUTBASE 'blue_'];
outbaseColor = [OUTBASE 'color_'];

% Generate output directory names
outdirRed = fullfile(OUTDIR, 'red');
outdirGreen = fullfile(OUTDIR, 'green');
outdirBlue = fullfile(OUTDIR, 'blue');
outdirColor = fullfile(OUTDIR, 'color');

%%%%%%%%%%%%%%%%
% % CREATE DIRECTORIES %%
%%%%%%%%%%%%%%%%

% Create output parent directory if it doesn't exist already
if ~exist(OUTDIR, 'dir') % Check existence of directory
    fprintf(1, ['Created output parent directory at ' OUTDIR '\n']); % Inform the user
    mkdir(OUTDIR); % Create the directory
end

% Make red output director if it doesn't already exist
if ~exist(outdirRed, 'dir')% Check existence of directory
    fprintf(1, ['Created red-channel output directory at ' outdirRed '\n']);% Inform the user
    mkdir(outdirRed); % Create the directory
end

if nChannels > 1 % Create color directories if the images are color images
% Make green output director if it doesn't already exist
    if ~exist(outdirGreen, 'dir')% Check existence of directory
        fprintf(1, ['Created green-channel output directory at ' outdirGreen '\n']); % Inform the user
        mkdir(outdirGreen); % Create the directory
    end

% Make blue output director if it doesn't already exist
    if ~exist(outdirBlue, 'dir')% Check existence of directory
        fprintf(1, ['Created blue-channel output directory at ' outdirBlue '\n']); % Inform the user
        mkdir(outdirBlue); % Create the directory
    end

% Make color output director if it doesn't already exist
    if ~exist(outdirColor, 'dir') % Check existence of directory
        fprintf(1, ['Created color output directory at ' outdirColor '\n']); % Inform the user
        mkdir(outdirColor); % Create the directory
    end
end

% Force the kernel radius to be even, so that the kernel length is odd.
% This ensures that the kernel is symmetric (or antisymmetric) about its anchor point.
% If an odd or non-integer radius is input, the kernel radius is forced to the next
% largest even integer.
kernelRadius = 2 * ceil(KERNELRADIUS / 2); 

% Force kernelRadius = 4 if parallel processing is selected, because
% currently the parallel doesn't work with arbitrary kernel lengths. This
% is a bad bug and needs to be fixed.
if NPROCESSORS > 1
    kernelRadius = 4;
end

% Make raw image filepaths
mirrorPathsTop = flipud(paths(2:kernelRadius + 1, :));
mirrorPathsBottom = flipud(paths(end - kernelRadius: end - 1, :));
mirroredPaths = cat(1, mirrorPathsTop, paths, mirrorPathsBottom); % Concatonate paths to include reflected images

% Calculate the size of sliding image stack
stackSize = 2 * kernelRadius + 1;

% Create the initial stack of images for hessian transformation (red channel)
slidingStackRed = loadImages(mirroredPaths(1 : stackSize, :), 1);

% Select the initial anchor image
anchorImageRed = slidingStackRed(:, :, kernelRadius + 1); 

% Only deal with other channels if images are multi-channel images
if nChannels > 1
    slidingStackGreen = loadImages(mirroredPaths(1: stackSize, :), 2);
    slidingStackBlue = loadImages(mirroredPaths(1 : stackSize, :), 3);
    anchorImageGreen = slidingStackGreen(:, :, kernelRadius + 1); % Select the initial anchor image
    anchorImageBlue = slidingStackBlue(:, :, kernelRadius + 1); % Select the initial anchor image
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INITIALIZE PARALLEL PROCESSING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize a parallel processing pool if more than one processor is specified
if NPROCESSORS > 1
    if matlabpool('size') > 0; % Close any open parallel processing pools
        matlabpool close 
    end
   matlabpool(NPROCESSORS); % Start a new parallel processing pool
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PERFORM HESSIAN TRANSFORMATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% outBaseRed = [outbaseRed 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form output basename
% outpathsRed = filepaths(outdirRed, outBaseRed, NDIGITS, EXT, START, STOP); % Determine paths to saved images
% 
% outBaseGreen = [outbaseGreen 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form output basename
% outpathsGreen = filepaths(outdirGreen, outBaseGreen, NDIGITS, EXT, START, STOP); % Determine paths to saved images
% 
% outBaseBlue = [outbaseBlue 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form output basename
% outpathsBlue = filepaths(outdirBlue, outBaseBlue, NDIGITS, EXT, START, STOP); % Determine paths to saved images

% Distinguish between single-processor and parellel processing jobs
if NPROCESSORS > 1 
% Calculate hessian eigenimages of all specified images (parallel processing)
    parfor p = 1:nImages
        imageNumber = num2str(START + p - 1, ['%0' num2str(NDIGITS) '.0f']);
        
% Inform the user
        fprintf(1, 'Calculating hessian eigen-image for %s\n', [IMBASE num2str(imageNumber) EXT]);
        
        redAnchor = imreadchannel(paths(p, :), 1);

% Create the sliding stack of  (red channel) (need a better way to do this,
% but can't figure out the parfor restrictions)
        redStack = [

            imreadchannel(mirroredPaths(p, :), 1); ...
            imreadchannel(mirroredPaths(p + 1, :), 1); ...
            imreadchannel(mirroredPaths(p + 2, :), 1); ...
            imreadchannel(mirroredPaths(p + 3, :), 1); ...
            imreadchannel(mirroredPaths(p + 4, :), 1); ...
            imreadchannel(mirroredPaths(p + 5, :), 1); ...
            imreadchannel(mirroredPaths(p + 6, :), 1); ...
            imreadchannel(mirroredPaths(p + 7, :), 1); ...
            imreadchannel(mirroredPaths(p + 8, :), 1); ...

            ];

%   Perform the hessian transformation. The permute(reshape(...)) mess is
%   due to a restriction in parfor loops against passing around matrices
%   that are larger than [m x n x 1]. The image stack above comes out as a
%   large [m x n], but hessian3 expects a 3D matrix for "stack," so all of
%   the reshaping has to be done jwithin the function call.
        hessian3(redAnchor(:, :, 1), permute(reshape(redStack', width, height, stackSize), [2 1 3]), GAUSSIANSCALES, RESCALE, SCALINGCOEFFICIENTS, outdirRed, outbaseRed, imageNumber, EXT);
        
% Only deal with the other channels if they exist
        if nChannels > 1           

% Select green anchor image
        greenAnchor = imreadchannel(paths(p, :), 1);

% Create the sliding stack of  (red channel) (need a better way to do this,
% but can't figure out the parfor restrictions)
        greenStack = [

            imreadchannel(mirroredPaths(p, :), 2); ...
            imreadchannel(mirroredPaths(p + 1, :), 2); ...
            imreadchannel(mirroredPaths(p + 2, :), 2); ...
            imreadchannel(mirroredPaths(p + 3, :), 2); ...
            imreadchannel(mirroredPaths(p + 4, :), 2); ...
            imreadchannel(mirroredPaths(p + 5, :), 2); ...
            imreadchannel(mirroredPaths(p + 6, :), 2); ...
            imreadchannel(mirroredPaths(p + 7, :), 2); ...
            imreadchannel(mirroredPaths(p + 8, :), 2); ...

            ];

%   Perform the hessian transformation. The permute(reshape(...)) mess is
%   due to a restriction in parfor loops against passing around matrices
%   that are larger than [m x n x 1]. The image stack above comes out as a
%   large [m x n], but hessian3 expects a 3D matrix for "stack," so all of
%   the reshaping has to be done jwithin the function call.
            hessian3(greenAnchor(:, :, 1), permute(reshape(greenStack', width, height, stackSize), [2 1 3]), GAUSSIANSCALES, RESCALE, SCALINGCOEFFICIENTS, outdirGreen, outbaseGreen, imageNumber, EXT);

% Select green anchor image
        blueAnchor = imreadchannel(paths(p, :), 1);

% Create the sliding stack of  (red channel) (need a better way to do this,
% but can't figure out the parfor restrictions)
        blueStack = [

            imreadchannel(mirroredPaths(p, :), 3); ...
            imreadchannel(mirroredPaths(p + 1, :), 3); ...
            imreadchannel(mirroredPaths(p + 2, :), 3); ...
            imreadchannel(mirroredPaths(p + 3, :), 3); ...
            imreadchannel(mirroredPaths(p + 4, :), 3); ...
            imreadchannel(mirroredPaths(p + 5, :), 3); ...
            imreadchannel(mirroredPaths(p + 6, :), 3); ...
            imreadchannel(mirroredPaths(p + 7, :), 3); ...
            imreadchannel(mirroredPaths(p + 8, :), 3); ...

            ];

%   Perform the hessian transformation. The permute(reshape(...)) mess is
%   due to a restriction in parfor loops against passing around matrices
%   that are larger than [m x n x 1]. The image stack above comes out as a
%   large [m x n], but hessian3 expects a 3D matrix for "stack," so all of
%   the reshaping has to be done jwithin the function call.
            hessian3(blueAnchor(:, :, 1), permute(reshape(blueStack', width, height, stackSize), [2 1 3]), GAUSSIANSCALES, RESCALE, SCALINGCOEFFICIENTS, outdirBlue, outbaseBlue, imageNumber, EXT);   

% Combine color channels
%             fprintf(1, '\n\nCombining color channels ... \n'); % Inform the user
            currentImage = START + p - 1; % Number of the current working image
% Combine channels for rescaled cases
            if RESCALE
                for h = 1:length(SCALINGCOEFFICIENTS)
                    outBaseRed = [outbaseRed 'scale' num2str(SCALINGCOEFFICIENTS(h)) '_'];
                    outBaseGreen = [outbaseGreen 'scale' num2str(SCALINGCOEFFICIENTS(h)) '_'];
                    outBaseBlue = [outbaseBlue 'scale' num2str(SCALINGCOEFFICIENTS(h)) '_'];
                    outBaseColor = [outbaseColor 'scale' num2str(SCALINGCOEFFICIENTS(h)) '_'];
                    combineChannels(outdirRed, outdirGreen, outdirBlue, outBaseRed, outBaseGreen, outBaseBlue, outdirColor, outBaseColor, NDIGITS, EXT, currentImage, currentImage);
                end
            else

% Combine channels for no rescaling
                combineChannels(outdirRed, outdirGreen, outdirBlue, outbaseRed, outbaseGreen, outbaseBlue, outdirColor, outbaseColor, NDIGITS, EXT, currentImage, currentImage);
            end
        end

    end

else

% Calculate hessian eigenimages of all specified images (single processor)
    for p = 1:nImages
        imageNumber = num2str(START + p - 1, ['%0' num2str(NDIGITS) '.0f']); % Current image number
        
        % Inform the user
        fprintf(1, 'Calculating hessian eigen-image for %s\n', [IMBASE num2str(imageNumber) EXT]);
                
        % Calculate and write hessian image
        hessian3(anchorImageRed, slidingStackRed, GAUSSIANSCALES, RESCALE, SCALINGCOEFFICIENTS, outdirRed, outbaseRed, imageNumber, EXT);

        if p < nImages % Don't advance the counter on the last step      
% Slide the image stack to prepare for the next hessian transformation 
            lastImageRed = loadImages(mirroredPaths(p + stackSize, :), 1); % Load the new anchor image
            slidingStackRed = circshift(slidingStackRed, [0 0 -1]); % Shift the elements of the image stack
            slidingStackRed(:, :, stackSize) = lastImageRed; % Append new image to sliding stack (red)
            anchorImageRed = slidingStackRed(:, :, kernelRadius + 1);  % Select the new anchor image (red)            
        end

% Only deal with the other channels if they exist
        if nChannels > 1

% Hessian transform of green channel
            hessian3(anchorImageGreen, slidingStackGreen, GAUSSIANSCALES, RESCALE, SCALINGCOEFFICIENTS, outdirGreen, outbaseGreen, imageNumber, EXT);

% Hessian transform of blue channel
            hessian3(anchorImageBlue, slidingStackBlue, GAUSSIANSCALES, RESCALE, SCALINGCOEFFICIENTS, outdirBlue, outbaseBlue, imageNumber, EXT);

            if p < nImages % Don't advance the counter on the last step
                lastImageGreen = loadImages(mirroredPaths(p + stackSize, :), 2); % Load the new anchor image (green)
                lastImageBlue = loadImages(mirroredPaths(p + stackSize, :), 3); % Load the new anchor image (blue)

                slidingStackGreen = circshift(slidingStackGreen, [0 0 -1]); % Shift the elements of the image stack (green)
                slidingStackBlue = circshift(slidingStackBlue, [0 0 -1]); % Shift the elements of the image stack (blue)

                slidingStackGreen(:, :, stackSize) = lastImageGreen; % Append new image to sliding stack (green)
                slidingStackBlue(:, :, stackSize) = lastImageBlue; % Append new image to sliding stack (blue)

                anchorImageGreen = slidingStackGreen(:, :, kernelRadius + 1);  % Select the new anchor image (green)
                anchorImageBlue = slidingStackBlue(:, :, kernelRadius + 1);  % Select the new anchor image (blue)
            end

% % Combine color channels
%             fprintf(1, '\n\nCombining color channels ... \n'); % Inform the user
            currentImage = START + p - 1; % Number of the working image
% Combine channels for rescaled cases
            if RESCALE
                for h = 1:length(SCALINGCOEFFICIENTS)
                    outBaseRed = [outbaseRed 'scale' num2str(SCALINGCOEFFICIENTS(h)) '_']; % Red channel output basename
                    outBaseGreen = [outbaseGreen 'scale' num2str(SCALINGCOEFFICIENTS(h)) '_']; % Green channel output basename
                    outBaseBlue = [outbaseBlue 'scale' num2str(SCALINGCOEFFICIENTS(h)) '_']; % Blue channel output basename
                    outBaseColor = [outbaseColor 'scale' num2str(SCALINGCOEFFICIENTS(h)) '_']; % Color channel output basename

                    % Combine the channels and write the output
                    combineChannels(outdirRed, outdirGreen, outdirBlue, outBaseRed, outBaseGreen, outBaseBlue, outdirColor, outBaseColor, NDIGITS, EXT, currentImage, currentImage);
                end
            else
                
% Combine channels for no rescaling
                combineChannels(outdirRed, outdirGreen, outdirBlue, outbaseRed, outbaseGreen, outbaseBlue, outdirColor, outbaseColor, NDIGITS, EXT, currentImage, currentImage);
            end
            
        end

    end

end
    
% Close any open parallel processing pools
if matlabpool('size') > 0
    matlabpool close; % Close parallel processing session
end

end

%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%%




