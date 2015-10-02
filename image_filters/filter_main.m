%% Image Processing Filter Suite

%% Available Filters:
% background_subtraction: Subtracts the temporal mean, median or minimum of the 
%   entire image set from each individual image. Because the filter is
%   performed on integers, no scaling is applied.
% 
% temporal_bandpass: Independently filters each pixel signal over time with
%   windowed Gaussian filters in the frequency domain. Because the filter
%   requires and outputs images of class double, a global scaling function
%   is applied after filtering to convert the data to 16-bit integer images 
%   scaled to the the maximum and minimum of the entire data set.
% 
% spatial_highpass: Independently filters each image in both dimensions with
%   windowed Gaussian filters in the frequency domain. Because the filter
%   requires and outputs images of class double, a global scaling function
%   is applied after filtering to convert the data to 16-bit integer images 
%   scaled to the the maximum and minimum of the entire data set.
% 
% crop_and_rotate: Rotates the image using Matlab imrotate to align flow
%   with Cartesian coordinates and crops using imcrop to the desired ROI.
%   Because the filter is performed on integers, no scaling is applied.
% 
% intensity_threshold: Thresholds the data intensity range to enhance
%   contrast of middle values. Only supports images with class uint8 or
%   uint16.
%
% phase_retrieval: Performs phase retrieval as described by Fouras group
%   (Paganin 2004, Irvine 2008). Because the filter requires and outputs 
%   images of class double, a global scaling function is applied after
%   filtering to convert the data to 16-bit integer images scaled to the 
%   the maximum and minimum of the entire data set.
%
% submedian_inversion: Subtracts the median value of the image from each
%   pixel and computes the absolute value of the image. Designed to convert
%   dark diffraction rings to positive signal. Works best after background
%   subtraction.
% highpass_butterworth: Transforms the image set to the frequency domain
%   and performs a frequency cutoff of the intensities.
%
% This file was editted on 12/8/2011 by Brett Meyers


%% Filter and directory setup
% A. 2 filter structure options are available:
%   1. additive = 1 - Filters will be applied sequentially, each filter to 
%   the images in the previously written directory.
%   2. additive = 0 - Filters will be applied separately and only to the
%   images in the parent directory.
Data.additive           = 1;
% B. Filters will be applied in the order listed. Separate strings by
% semicolons. Available filters are listed in the previous cell.
% Data.filter_to_run = {'phase_retrieval';'submedian_inversion'};
Data.filter_to_run = {'background_subtraction'};
% C. 3 write directory options are available:
%   1. write_directory = [] - A cascading file structure will be created in
%   the parent directory. Images output by each filter will be written in a
%   subfolder of the previous directory.
%   2. write_directory = 'name' - This is similar to the first option but
%   instead of creating the file cascade in the parent directory, it will be
%   created in the directory 'name'. End directory name with \ or /.
%   3. write_directory = {'name1';'name2';...;'nameN'} - This will write the
%   output of each filter to the corresponding write_directory named by the
%   user. The number of entries in write_directory must be equal to the
%   number of filters called (number of entries in filter_to_run). End
%   directory names with \ or /.
Data.write_directory    = [];

%% Filter Parameters
% Parameters required for all filters
% Data.parent_directory   = '/Volumes/current_storage/Projects/WFU/echoPIV/patients/PIV-066/TIF/01--02-42-00-4c/'; % Directory with original images, must end with \ or /5
% Data.parent_directory   = '/Users/brettmeyers/Google Drive/test_images/nh_test_01/red/'; % Directory with original images, must end with \ or /5
Data.parent_directory = '/Users/brettmeyers/Desktop/Mouse_Data/images/Long Axis1 Sharp/';
Data.image_format       = 'tif';
Data.frame_min          = 1;     % index of first frame in directory for computation
Data.frame_max          = Inf;    % index of last frame in directory for computation, set Inf for last frame in directory

% Filter-specific parameters - see filter functions for further details
% background_subtraction
Data.background_subtraction.to_subtract = 'mean';               % 'mean' or 'median' or 'minimum

% temporal_bandpass
Data.temporal_bandpass.sigma_vect = [0.02 0.08];                % Gaussian standard deviation for lowpass and highpass filter respectively, between 0 and 1 (unitless)

% spatial_highpass
Data.spatial_highpass.sigma = 0.025;                             % Gaussian standard deviation for highpass filter, between 0 and 1

% crop_and_rotate
Data.crop_and_rotate.angle = -90.5;                               % angle of rotation (degrees)
Data.crop_and_rotate.crop_window = [220 400 670 800];            % [x0 y0 width height] (pixels)

% intensity_threshold
Data.intensity_threshold.threshold = [0.1 0.3];              % lower and upper intensity thresholds respectively, between 0 and 1 (fraction of maximum intensity of image set)

% phase_retrieval
Data.phase_retrieval.tau = 10;                                 % filter tuning variable (unitless)

%highpass_butterworth
Data.n_order = 6;
Data.cutoff_freq = 0.1;
Data.ftype = 'high';

%median window subtraction
Data.median_window_sub.to_subtract = 'mean';               % 'mean' or 'median' or 'minimum
Data.windowsize=3;

%mask region removal
Data.windowdim=3;

%discrete wavelet transform
Data.levelnumber = 2;
Data.waveletname = 'sym8';
Data.alpha       = 3;
Data.gbl_or_lvd  = 'lvd';
Data.s_or_h      = 's';
Data.keep_approx = 1;


%% Auto generation of write directories
Data.write_directory = generate_write_directories(Data);

%% Auto filter call
fprintf('-- Image Processing --\n\n')
Data.read_directory = Data.parent_directory;
% iterate through filters
for i=1:size(Data.filter_to_run)
    Data.index = i;
    fprintf('Running filter %s\n',Data.filter_to_run{i})
    
    eval([Data.filter_to_run{i},'(Data);'])                                 % run filter
    
    save(fullfile(Data.write_directory{i},'filter_main_Data.mat'),'Data')   % save variables to output directory
    
    if Data.additive                                                        % update read directory if filters are additive
        Data.read_directory = Data.write_directory{i};
    end
    fprintf('\n')
end
fprintf('\n')