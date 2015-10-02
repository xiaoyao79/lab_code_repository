% Use this template for batch-processing hessian jobs.

IMDIR = fullfile(pwd, '..', 'testImages');  % Path to directory containing images
IMBASE = 'testImage_'; % Image base name (i.e., IMBASE = 'image_' for image_00001.tif)
OUTDIR = '~/Desktop/hessian_test'; % Path to directory UNDERNEATH WHICH the output directories will be created.
OUTBASE = 'testImage_hessian_'; % Base name of output images
NDIGITS = 3; % Number of digits in the file name (i.e., NDIGITS = 5 for image_00001.tif)
EXT = '.tif'; % Image extension (i.e., EXT = '.tif' for image_00001.tif)
START = 1; % Number of first image to be processed
STOP = 20; % Number of final image to be processed

GAUSSIANSCALES = [1 1 1]; % Standard deviations (in pixels) in the [ X Y Z ] directions of the gaussian kernel used for differentiation
KERNELRADIUS = 4; % Radius of the kernel used for the differentiation scheme. This is forced to 4 for parallel processing for now.
RESCALE = 0; % Binary flag that specifies whether or not to rescale the magnitudes of the time derivaties.
NPROCESSORS = 7; % Number of processors to use. 
SCALINGCOEFFICIENTS = 10; % Target ratio of rescaled time-derivative to X and Y derivatives. Only used when RESCALE = 1. 
% If SCALINGCOEFFICIENTS is a vector, then hessianRUN will perform the
% hessian operation with scales corresponding to each element of
% SCALINGCOEFFICIENTS.

% Run the Hessian code
hessianRun(IMDIR, IMBASE, OUTDIR, OUTBASE, NDIGITS, EXT, START, STOP, GAUSSIANSCALES, KERNELRADIUS, RESCALE, SCALINGCOEFFICIENTS, NPROCESSORS)


