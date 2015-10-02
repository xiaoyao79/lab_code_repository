function[Xf] = denoise_echopiv(IM1)
%% Discrete Wavelet Transform for De-noising and Sharpening EchoPIV Images
%  Current code takes spatial (2-D) inputs for each image and performs
%  a series of operations to de-noise and sharpen images.  This will aide
%  in improved image correlation during processing, resulting in more
%  accurate flowfields.  
%  
%  Written by Brett Meyers on 2/3/2012
%
% Required input parameters for wavelet denoising:
%   IM:     Image going through denoising process
%   N:      Wavelet decomposition level
%   wname:  Mother wavelet being selected for decomposition
%   thr:    Selected threshold for de-noising.
%   de_comp:Selection of denoising, 'gbl', or compresion, 'lvd'
%   SORH:   Soft or hard thresholding
%   KEEPAPP:Perfoms to either threshold or preserve approximation
%       coefficients. If KEEPAPP = 1 no thresholding of app. coef. will
%       occur.
%
% Returned outputs from algorithm:
%   c:      Decomposition vector of image
%   s:      Bookkeeping matrix of decomposition
%   A:      Matrix of approximation co-efficients required to find
%           threshold levels
%   thr:    Selected threshold for de-noising
%   xc:     Denoised version of image
%   perf0:  L^2 norm recovery score
%   perf2:  L^2 norm compression score

% Input parameters
IM      = im;
N       = 2;
wname   = 'sym8';
alpha   = 3;
%thr     = 20;
de_comp = 'lvd';
SORH    = 'h';
KEEPAPP = 1;

% Performs wavelet decomposition on image, calculates the approximate
% coeffictient, and returns an image threshold for de-noising.
[c,s] = wavedec2(IM,N,wname);
A = appcoef2(c,s,wname,1);
%thr = thselect(A,'rigrsure');
[thr,nkeep] = wdcbm2(c,s,alpha,6*prod(s(1,:)));

% Performs de-noising and reconstructs the now de-noised image from the
% wavelet information.
[xc,cxc,lxc,perf0,perf12] = wdencmp(de_comp,c,s,wname,N,thr,SORH);
X = uint8(xc);

% Displays filtered image versus unfiltered image
figure(1)
imshow(X);
figure(2)
imshow(IM);
figure(3)
% keyboard;
imshow(imadjust(IM-X));

w  = fspecial('unsharp');
Xs = single(X);
%w  = [[1,1,1];[1,-8,1];[1,1,1]];
Xf = uint8(imfilter(Xs,w,'replicate'));
figure(4)
imshow(Xf)
figure(5)
imshow(X-Xf)
