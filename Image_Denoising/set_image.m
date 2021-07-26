function [u_exact, u0] = set_image(noiseLevel,DimensionSize)
% the script is to read the camera-man image and add noise to it with
% the corresponding Gaussian mean and variance. 
image_exact = imread('cameraman.tif');
% re-scale the image to dimensionSize of problem and convert to double
% precision: 
u_exact = im2double(imresize(image_exact, DimensionSize)); 


% add random noise to the image: 
randn('seed',0);
u0 = u_exact+noiseLevel*randn(size(u_exact));
u0 = (u0-min(min(u0)))/(max(max(u0))-min(min(u0)));
end 