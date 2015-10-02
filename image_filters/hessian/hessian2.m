function IMAGE = hessian2(IMAGE, SCALES)
% hessian2(IMAGE) Transforms an IMAGE such that the pixel values in the new image are the
% greatest magnitudes of the eigenvalues of the 2 x 2 intensity gradient matrix (the Hessian matrix).
% 
% hessian2(IMAGE, SCALES) performs the same transformation, with the
% standard deviation of the gaussian convolution kernel used calculate derivatives
% multiplied in each direction by the elements of the vector SCALE, where SCALE = [ scaleX scaleY ].
% The default SCALES vector is [ 1 1 ], i.e., the standard deviation of the
% convolution kernel is 1 pixel in each direction.
%
% KNOWN ISSUES
%    If you're running this function as part of a loop, you may increase speed
%   by disabling the print-to-screen status updates. This isn't implemented as a
%   binary flag anywhere right now, but for the time being, just comment out
%   the "fprintf(1, ...)" commands as you see fit. 
% 
% INPUTS
%   IMAGE = 2-D array containing image intensity data. The [ 1, 2 ] dimensions of the array
%   correspond to the [ y, x ] dimensions of the image image.
%
%   SCALES = 2 - element vector containing the standard deviations of the
%   gaussian kernels used to calculate image derivatives in each direction.
%   The elements of SCALE are [ scaleX scaleY ]. The default value
%   of SCALES is [ 1 1 ] .
% 
% OUTPUTS
%   IMAGE = Image whose pixel values are  the values of the
%   greatest-magnitude eigenvalues of the 2 x 2 intensity gradient matrix
%   (the Hessian matrix). The variable IMAGE is reused here to save memory.
%   
% SEE ALSO
%   imdifferentiate

%%%%%%%%%%%%%%%
%%% BEGIN FUNCTION %%%
%%%%%%%%%%%%%%%

% Default to scales of 1 in every direction.
% SCALES = [ScaleX ScaleY]
if nargin < 2
    SCALES = [ 1 1 ]; 
end

% Throw an error if image is not a 2-D image
if ndims(IMAGE) ~= 2
    error('Images must be 2-D images');
end

% Calculate the dimensions of the image
[height width] = size(IMAGE);

% Second derivatives (diagonal terms)
fprintf(1, 'Calculating Ixx ... \n');
dxx = single(imdifferentiate(IMAGE, SCALES, 2, 0, 0)); % Second derivative with respect to X twice. Store in single-precision format.
fprintf(1, 'Calculating Iyy ... \n');
dyy = single(imdifferentiate(IMAGE, SCALES, 0, 2, 0)); % Second derivative with respect to Y twice Store in single-precision format.

% Second derivatives (cross term)
fprintf(1, 'Calculating Ixy ... \n');
% Second derivative with respect to X then Y. Store in single-precision format.
dxy = single(imdifferentiate(IMAGE, SCALES, 1, 1, 0 )); 

% Inform the user
fprintf(1, 'Calculating coefficients of the characteristic equation ... \n');

% Calculate the coefficients of the eigenvalue characteristic equation L^2 + b*L + c = 0
% Coefficients on the L term 
b = - (dxx + dyy);

% Constants not multiplied by L
c = (dxx.*dyy) - (dxy.*dxy);

% Sign of the coefficient b in the characteristic equation
signs = -1 * ones(height, width) + 2 * (b > 0); 

% Inform the user
fprintf(1, 'Calculating roots of the characteristic equation ... \n');

% Roots of the characteristic equation
q = - 0.5 * (b + signs .* sqrt(b .^2 - 4 *c)); 

% Maximum eigenvalue, as calculated by ImageJ. Not sure why c is divided by
% q in the ImageJ code, but the results look reasonable. This may change in
% future revisions of this MATLAB function.
IMAGE = max(abs(q), abs(c./q)); 

% Inform the user
fprintf(1, 'Done! \n');

end

%%%%%%%%%%%%%%%%
%%% END OF FUNCTION % % %
%%%%%%%%%%%%%%%%






