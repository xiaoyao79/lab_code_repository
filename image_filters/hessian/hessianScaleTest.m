function hessianScaleTest(IMDIR, IMBASE, OUTDIR, OUTBASE, NDIGITS, EXT, START, STOP, SCALINGCOEFFICIENTS, NPROCESSORS)
% hessianScaleTest(IMDIR, IMBASE, OUTDIR, OUTBASE, NDIGITS, EXT, START, STOP, SCALINGCOEFFICIENTS, NPROCESSORS) 
% Applies a range of z-derivative scaling coefficients to an image sequence
% and saves the results.
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
%   SCALINGCOEFFICIENTS = Vector containing values of forced ratios of
%   the mean magnitude of the Z- intensity derivative to those of the
%   X and Y intensity derivatives. 
%
%   NPROCESSORS = Number of processors to use in parallel (integer)
% 
% OUTPUTS
%   None (images are saved to disk)
%
% SEE ALSO
%   hessian3, checkHessianJobMemory, filepaths, combineChannels

%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%


% Default to 1 processor
if nargin < 10
    NPROCESSORS = 1;
end

% These variables should be fixed for this test
AUTOSCALE = 1;
SCALES = [1 1 1];

% Determine whether job will run with the given system memory
% If the job memory exceeds 75% of the system memory, decrease the number
% of processors used. If the single-processor job still exceeds 75% of the
% system memory, decrease the number of images on which to perform the
% hessian transformation. 
[STOP NPROCESSORS] = checkHessianJobMemory(IMDIR, IMBASE, NDIGITS, EXT, START, STOP, NPROCESSORS);

% Make raw image filepaths
paths = filepaths(IMDIR, IMBASE, NDIGITS, EXT, START, STOP);

% Determine number of color channels
nChannels = size(imread(paths(1, :)), 3); 

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

% Create output directory if it doesn't exist already
if ~exist(OUTDIR, 'dir')
    fprintf(1, ['Created output parent directory at ' OUTDIR '\n']);
    mkdir(OUTDIR);
end

% Make red output director if it doesn't already exist
if ~exist(outdirRed, 'dir')
    fprintf(1, ['Created red-channel output directory at ' outdirRed '\n']);
    mkdir(outdirRed);
end

if nChannels > 1 % if the images are color images
% Make green output director if it doesn't already exist
    if ~exist(outdirGreen, 'dir')
        fprintf(1, ['Created green-channel output directory at ' outdirGreen '\n']);
        mkdir(outdirGreen);
    end

% Make blue output director if it doesn't already exist
    if ~exist(outdirBlue, 'dir')
        fprintf(1, ['Created blue-channel output directory at ' outdirBlue '\n']);
        mkdir(outdirBlue);
    end

% Make color output director if it doesn't already exist
    if ~exist(outdirColor, 'dir')
        fprintf(1, ['Created color output directory at ' outdirColor '\n']);
        mkdir(outdirColor);
    end
end

% Close any open parallel processing pools
if matlabpool('size') > 0; 
    matlabpool close
end

% Start a new parallel processing pool
matlabpool(NPROCESSORS)

% Compute hessian eigenimages
% Red channels
fprintf(1, '\n\nCalculating eigenimages for RED channels ... \n'); % Inform the user
fprintf(1, 'Loading raw images... \n'); % Inform the user
rawImages = loadImages(paths, 1); % Load red channel of raw images
parfor k = 1:length(SCALINGCOEFFICIENTS)
    fprintf(1, ['Calculating eigenimages for SCALE = ' num2str(SCALINGCOEFFICIENTS(k)) '\n']); % Inform the user
    outBase = [outbaseRed 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form output basename
    outpaths = filepaths(outdirRed, outBase, NDIGITS, EXT, START, STOP); % Determine paths to output images
    IMAGE = hessian3(rawImages, SCALES, AUTOSCALE, SCALINGCOEFFICIENTS(k)); % Calculate hessian eigenimages
    for m = 1:size(outpaths, 1); imwrite(uint8(IMAGE(:, :, m)), outpaths(m, :)); end; % Write images to disk
end

if nChannels > 1 % If the images are color images

% Green channels
    fprintf(1, '\n\nCalculating eigenimages for GREEN channels ... \n'); % Inform the user
    fprintf(1, 'Loading raw images... \n'); % Inform the user
    rawImages = loadImages(paths, 2); % Load green channel of raw images
    parfor k = 1:length(SCALINGCOEFFICIENTS)
        fprintf(1, ['Calculating eigenimages for SCALE = ' num2str(SCALINGCOEFFICIENTS(k)) '\n']); % Inform the user
        outBase = [outbaseGreen 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form output basename
        outpaths = filepaths(outdirGreen, outBase, NDIGITS, EXT, START, STOP); % Determine paths to output images
        IMAGE = hessian3(rawImages, SCALES, AUTOSCALE, SCALINGCOEFFICIENTS(k)); % Calculate hessian eigenimages
        for m = 1:size(outpaths, 1); imwrite(uint8(IMAGE(:, :, m)), outpaths(m, :)); end; % Write images
    end

% Blue channels
    fprintf(1, '\n\nCalculating eigenimages for BLUE channels ... \n'); % Inform the user
    fprintf(1, 'Loading raw images... \n'); % Inform the user
    rawImages = loadImages(paths, 3); % Load blue channel of raw images
    parfor k = 1:length(SCALINGCOEFFICIENTS)
        fprintf(1, ['Calculating eigenimages for SCALE = ' num2str(SCALINGCOEFFICIENTS(k)) '\n']); % Inform the user
        outBase = [outbaseBlue 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form output basename
        outpaths = filepaths(outdirBlue, outBase, NDIGITS, EXT, START, STOP); % Determine paths to output images
        IMAGE = hessian3(rawImages, SCALES, AUTOSCALE, SCALINGCOEFFICIENTS(k)); % Calculate hessian eigenimages
        for m = 1:size(outpaths, 1); imwrite(uint8(IMAGE(:, :, m)), outpaths(m, :)); end; % Write images  
    end

% Combine color channels
    fprintf(1, '\n\nCombining color channels ... \n'); % Inform the user
    parfor k = 1:length(SCALINGCOEFFICIENTS)
        fprintf(1, ['Combining color channels for SCALE = ' num2str(SCALINGCOEFFICIENTS(k)) '\n']); % Inform the user
        imbaseRed = [outbaseRed 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form red channel basename
        imbaseGreen = [outbaseGreen 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form red channel basename
        imbaseBlue = [outbaseBlue  'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form red channel basename
        outBase = [outbaseColor 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_']; % Form output basename
        combineChannels(outdirRed, outdirGreen, outdirBlue, imbaseRed, imbaseGreen, imbaseBlue, outdirColor, outBase, NDIGITS, EXT, START, STOP)
    end

end

matlabpool close; % Close parallel processing session

end

%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%%




