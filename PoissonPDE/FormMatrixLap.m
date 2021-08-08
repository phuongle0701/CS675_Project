% function to implement the Laplacian optimisation problem: 
function [A,b] = FormMatrixLap(u,alpha)
n = size(u,1);

Dn = eye(n); 
Dxx = diff(Dn,2,1); 
Dyy = diff(Dn,2,1); 
LA = kron(Dxx,Dn) + kron(Dn,Dyy); % laplacian operator. 


I = eye(n*n); 
A1 = zeros(n,n); 
A1(1:alpha,:) = 1; % alpha 
A1(n-alpha:n,:) = 1; 
B1 = zeros(n,n); % B is A' 
B1(1:alpha, 1) = 1; 
B1(n:n+alpha,:) = 1; 


A = [LA;sqrt(alpha)*I(A1,:)]; 
b = [zeros(size(LA,1),1); B(B1)]; 
end 