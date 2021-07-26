function driver 
clear all; 
clc; 
noiseLevel = 0.05; 
[exact, noise_img] = set_image(noiseLevel,[256 256]); 
set(gcf, 'Position',  [100, 100, 200, 300])
figure(1); 
montage({exact, noise_img}, 'Size', [1 2]); 
disp(size(exact)); 
end 