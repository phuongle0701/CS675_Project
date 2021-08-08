function check_lap
clear all; 
clc; 
format short; 

[A,b] = Laplacian(50,'2D'); 
x = A\b; 

X = reshape(x, 50,50); 
imagesc(X);
colorbar('AxisLocation','out');


end 