function IMAGE = loadImages(FILEPATHS, CHANNEL)
% loadImages(FILEPATHS, CHANNEL) loads a stack of images of channel CHANNEL.
%
% INPUTS
%   FILEPATHS = Matrix whose rows are strings containing the file paths to
%   each image in the specified sequence (strings). 
%
%   CHANNEL = Image channel to load (integer)
%
% OUTPUTS
%   IMAGE = [m -by- n -by- p] matrix containing images of height (number of
%   rows) m, width (number of columns)n, and number of images p. (uint8 integers)
%
% SEE ALSO
%   filepaths
%

%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%

% Default to red channel
if nargin < 2
    CHANNEL = 1;
end

nImages = size(FILEPATHS, 1); % Number of images
[height width nChannels] = size(imread(FILEPATHS(1, :)));% Read first image to determine its size

IMAGE = uint8(zeros(height, width, nImages)); % Initialize Image matrix


% If the specified channel exceeds the number of channels, load the final channel.
if CHANNEL > nChannels
    CHANNEL = nChannels;
end

% Load images
% h = waitbar(0, 'Loading images ... '); %Initialize wait bar
% fprintf(1, 'Loading images ... \n');
for k = 1:nImages
    tempImage = uint8(imread(FILEPATHS(k, :))); % Read image
    IMAGE(:, :, k) = tempImage(:, :, CHANNEL); % Populate IMAGE matrix with the specified channel of the image
%     waitbar(k/nImages); % Update wait bar
end

% close(h) % Close the waitbar figure

% pause(1/1000) % Pause briefly so that the waitbar figure closes

IMAGE = uint8(IMAGE); % Make sure IMAGE is uint8 format.

end

%%%%%%%%%%%%%%%%%%%%
%%% END FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%

