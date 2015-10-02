function IMAGE = hessian3(ANCHOR, STACK, GAUSSIANSCALES, RESCALE, SCALINGCOEFFICIENTS, SAVEDIR, SAVEBASE, IMAGENUMBER, EXT)
% IMAGE = hessian3(ANCHOR, STACK, GAUSSIANSCALES, RESCALE, SCALINGCOEFFICIENTS, SAVEDIR, SAVEBASE, IMAGENUMBER, EXT) Transforms an IMAGE such that the pixel values in the new image are the
% greatest magnitudes of the eigenvalues of the 3 x 3 intensity gradient matrix
% (the Hessian matrix).
% 
% hessian3(IMAGE, SCALES) performs the same transformation, with the
% standard deviation of the gaussian convolution kernel used calculate derivatives
% multiplied in each direction by the elements of the vector SCALE, where SCALE = [scaleX scaleY scaleZ].
% The default SCALES vector is [1 1 1], i.e., the standard deviation of the
% convolution kernel is 1 pixel in each direction.
%
% KNOWN ISSUES
%   Right now this code won't take 2-D images. I need to add an if statement
%   to calculate the eigenvalues of a  2 x 2 hessian matrix (i.e., for 2D images). 
% 
% INPUTS
%   ANCHOR = Anchor image the 3-D hessian operation. 
%
%   STACK = Stack of images for the 3-D hessian operation
%
%   GAUSSIANSCALES = 3 - element vector containing the standard deviations of the
%   gaussian kernels used to calculate image derivatives in each direction.
%   The elements of SCALE are [ scaleX scaleY scaleZ ]. The default value
%   of SCALES is [ 1 1 1 ] .
% 
%   RESCALE = Binary flag that specifies whether or not to scale the
%   Z-derivatives by a constant value such that the order of the magnitude of the  Z
%   derivatives is similar to those of the X- and Y- derivatives. 
%
%   SCALINGCOEFFICIENTS = Forced ratio of the mean magnitude of the
%   Z- intensity derivative to those of the X and Y intensity derivatives.
%   SCALINGCOEFFICIENT is used when RESCALE = 1. 
%
%   SAVEDIR = Path to directory in which output images are saved (string)
%
%   SAVEBASE = Basename of output images (string). For image_00001.tif,
%   SAVEBASE = 'image_'
%
%   IMAGENUMBER = Number of the image to be saved (formatted string). For image_00001.tif,
%   IMAGENUMBER = '00001';
%
%   EXT = File extension for the image to be saved (string). For image_00001.tif,
%   EXT = '.tif';
% 
% OUTPUTS
%   IMAGE = Image whose pixel values are  the values of the
%   greatest-magnitude eigenvalues of the 3 x 3 intensity gradient matrix
%   (the Hessian matrix). The variable IMAGE is reused here to save memory.
%   
% SEE ALSO
%   imdifferentiate

%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%
% Default to tiff images
if nargin < 9
    EXT = '.tif';
end

if nargin < 8;
    IMAGENUMBER = 0;
end

% Default to no saved image base name (i.e. don't save images)
if nargin < 7
    SAVEBASE = '';
end

% Default to no save path
if nargin < 6
    SAVEDIR = '';
end

% Default to scaling coefficient of 1;
if nargin < 5;
    SCALINGCOEFFICIENTS = 1;
end

% Default to not autoscale derivatives
if nargin < 4 
    RESCALE = 0;
end

% Default to scales of 1 in every direction.
% SCALES = [ScaleX ScaleY ScaleZ]
if nargin < 3
    GAUSSIANSCALES = [1 1 1]; 
end

% Size of image
dimensions = size(ANCHOR);

% Throw an error if image contains dimensions higher than 3
if length(dimensions) > 3
    error('Images must be less than or equal to 3 dimensions');
end

% Second derivatives (diagonal terms)
% fprintf(1, 'Calculating Ixx ... \n');
dxx = imdifferentiate(ANCHOR, GAUSSIANSCALES, 2, 0, 0); % Second derivative with respect to X twice. Store in single-precision format.
% fprintf(1, 'Calculating Iyy ... \n');
dyy = imdifferentiate(ANCHOR, GAUSSIANSCALES, 0, 2, 0); % Second derivative with respect to Y twice Store in single-precision format.
% fprintf(1, 'Calculating Izz ... \n');
dzz = imdifferentiate(STACK, GAUSSIANSCALES, 0, 0, 2); % Second derivative with respect to Z twice Store in single-precision format.

% Second derivatives (cross terms)
% fprintf(1, 'Calculating Ixy ... \n');
% Second derivative with respect to X then Y. Store in single-precision format.
dxy = imdifferentiate(ANCHOR, GAUSSIANSCALES, 1, 1, 0 ); 
% fprintf(1, 'Calculating Ixz ... \n');
% Second derivative with respect to X then Z. Scale according to the square root of the dzz scaling factor. Store in single-precision format.
dxz = imdifferentiate(STACK, GAUSSIANSCALES, 1, 0, 1); 
% fprintf(1, 'Calculating Iyz ... \n');
% Second derivative with respect to Y then Z. Scale according to the square root of the dzz scaling factor. Store in single-precision format.
dyz = imdifferentiate(STACK, GAUSSIANSCALES, 0, 1, 1); 

%Calculate scaling coefficient to force Z-derivatives to be of the same order as the X and Y derivatives
if RESCALE
    gradientRatio = 2 * mean(abs(dzz(:))) / ( mean(abs(dxx(:))) + mean(abs(dyy(:)))); % Calculate scaling coefficient
    scalingCoefficients = SCALINGCOEFFICIENTS / gradientRatio; % Scaling coefficient by which to multiply Z-derivatives
else
    scalingCoefficients = 1; % Set the scaling coefficient to 1 if RESCALE is not enabled.   
end

% Re-scale Z derivatives
for k = 1:length(SCALINGCOEFFICIENTS)
    dzzScaled = scalingCoefficients(k) * dzz; % Rescale Dzz derivative
    dxzScaled = sqrt(scalingCoefficients(k)) * dxz;% Rescale Dxz derivative
    dyzScaled = sqrt(scalingCoefficients(k)) * dyz; % Rescale Dyz derivative
    
% Compute the coefficients of the third-order polynomial whose roots are the
% eigenvalues the 3 x 3 intensity gradient matrix (the Hessian matrix).
% The polynomial is of the form L^3 + b * L^2 + c * L + d = 0, where L is
% the eivenvalue and b, c, and d are the coefficients derived from taking
% the determinaint of the 3 x 3 matrix in the eigenvalue problem.
% The eigenvalue characteristic equation used here is the actual 
% eivenvalue equation multiplied by -1 so that the coefficient
% on the L^3 term is equal to 1

% Inform the user
%     fprintf(1, 'Calculating coefficients of the characteristic equation ... \n');

% Coefficients on the L^2 term
    b = - ( dxx + dyy + dzzScaled );

% Coefficients on the L term
    c = dxx .* (dyy + dzzScaled) ...
     + dyy .* dzzScaled ...
     - dxy.*dxy ...
     - dyzScaled .* dyzScaled ...
     - dxzScaled .* dxzScaled;

% Constants in the eigenvalue polynomial
    d = dxx .* dyzScaled .* dyzScaled ...
      + dyy .* dxzScaled .* dxzScaled ...
      + dzzScaled .* dxy .* dxy ...
      - dxx .* dyy .* dzzScaled ...
      - 2 * (dxy .* dyzScaled .* dxzScaled);

% Substitute L = t - b/(3*a) the Tschirnhaus transformation, which puts the
% polynomial in the form t^3 + p*t + q = 0. Next we calculate the
% coeffieicnts p and q, with a = 1.

% Inform the user
%     fprintf(1, 'Calculating coefficients of the Tschirnhaus polynomial ... \n');

% Coefficient on the t term of the transformed polynomial
    p = (3 * c - b.^2) / 3; 

% Constant in the transformed polynomial
    q = (2 * b.^3 - 9 * b .* c + 27 * d.^2) / 27;

% Find the three roots trigonometrically. This method uses the substitution t = u * cos(theta) 
% to make the transformed polynomial coincide with the identity
% 4 *(cos(theta))^3  - 3 * cos(theta) - cos(3 * theta) = 0.
% Here u = 2 *sqrt( - p / 3). Real eigenvalues occur when p < 0.
%  The logical coefficient (1 - p > 0) sets the polynomial root to 0 if the root is imaginary. 

% %Inform the user
%     fprintf(1, 'Calculating first eigenvalue ... \n');

% Magnitute of first root
    t0 = real((1 - p > 0) .* sqrt((2 * sqrt(-p / 3) .* cos(1 / 3 * acos(3 * q ./ (2 * p) .* sqrt(- 3 ./ p)))).^2)); 

% %Inform the user
%     fprintf(1, 'Calculating second eigenvalue ... \n');

% Magnitute of second root
    t1 = real((1 - p > 0) .* sqrt((2 * sqrt(-p / 3) .* cos(1 / 3 * acos(3 * q ./ (2 * p) .* sqrt(- 3 ./ p)) - 2 * pi / 3)).^2));

% %Inform the user
%     fprintf(1, 'Calculating third eigenvalue ... \n');

% Magnitude of third root
    t2 = real((1 - p > 0) .* sqrt((2 * sqrt(-p / 3) .* cos(1 / 3 * acos(3 * q ./ (2 * p) .* sqrt(- 3 ./ p)) - 4 * pi / 3)).^2));

% % Inform the user
%     fprintf(1, 'Determining largest magnitude eigenvalues ... \n');

% Determine largest magnitude eigenvalues
    IMAGE = max(max(t0, t1), t2); 
    
% Deal with NaNs
    realImage = IMAGE(~isnan(IMAGE));
    
% Calculate the 95% intensity
    maxIntensity = mean(realImage(:)) + 2 * std(realImage(:));

% Normalize the hessian image to the 95% intensity
    IMAGE = IMAGE * 255 / maxIntensity;
    
% Convert the image to uint8 format
    IMAGE = uint8(IMAGE);
    
% Create path to saved file
    if RESCALE == 0; % If no rescaling was performed, don't include scaling information in the name of the saved image.
        savePath = fullfile(SAVEDIR, [SAVEBASE num2str(IMAGENUMBER) EXT]); 
    else
        
% If rescaling was performed, DO include scaling information in the name of
% the saved image
        savePath = fullfile(SAVEDIR, [SAVEBASE 'scale' num2str(SCALINGCOEFFICIENTS(k)) '_' num2str(IMAGENUMBER) EXT]);
    end
    
% Write images if save directory and base name were specified
    if (~strcmp(SAVEDIR, '') && ~strcmp(SAVEDIR, ''))  
        imwrite(IMAGE, savePath); % Writing image
        fprintf(1, 'Saved Hessian image to %s \n', savePath); % Inform the user
    end
    
% % Inform the user
%     fprintf(1, 'Done! \n');  
    
end

end

%%%%%%%%%%%%%%%%%%%%
%%% END FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%