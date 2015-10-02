function OUTPUTIMAGE = imreadchannel(PATH, CHANNEL)
% imreadchannel(PATH, CHANNEL) Reads an image located at PATH and returns a single color
% channel specified by CHANNEL
%
% INPUTS
%   PATH = File path to image (string)
%
%   CHANNEL = Channel of the image to return (integer, 1 ? CHANNEL ? 3)
% 
% OUTPUTS
%   OUTPUTIMAGE = Image channel exctracted from the image located at PATH (uint8)
%
% SEE ALSO
%   imread

%%%%%%%%%%%%%%%
%%% BEGIN FUNCTION %%%
%%%%%%%%%%%%%%%

% Read the image from disk
IMAGE = imread(PATH);

% Return the specified channel from the image
OUTPUTIMAGE = IMAGE(:, :, CHANNEL);

%%%%%%%%%%%%%%
%%% END FUNCTION%%%
%%%%%%%%%%%%%%

end