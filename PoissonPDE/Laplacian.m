% function to build the finite difference matrix for Laplacian PDE equation
% Input: 
% m: number of problem sizes. 
% dim: '1D', '2D', '3D' dimension of problem up to 3D dimension. 
% Output: matrix A represents the finite difference for Laplacian matrix. 
% b: matrix of the right hand side. 

% We use the approach of 4.8 Problem Van Loan-Golub 

% PDE equation: 
% -Delta(U) = f(X,Y,Z) 
% f(X,Y,Z): two central heating seating: 
% f = 1 if || ({x,y,z}) - ({0.3,0.3,0.3}) || <= 0.1 
% f = 1 if if || ({x,y,z}) - ({0.6,0.6,0.6}) || <= 0.1 
% f = 0 otherwise. 

function [A,b] = Laplacian(m, dim)
switch dim
    case '1D'
        [A,b] = Lap1D(m); 
    case '2D'
        [A,b] = Lap2D(m); 
    case '3D'
        [A,b] = Lap3D(m); 
end 
end 


% 1D laplacian
function [A,b] = Lap1D(N)
 
h = 1/(n+1); 
e = ones(n,1); 
% Build the stencil matrix 1D. 
A1D = spdiags([e -2*e e], -1:1, n,n);
% Build the left hand side matrix of finite difference.
A = -(1/(h^2))*A1D; 
% Store as sparse matrix: 
A = sparse(A); 
% build the heat-source: 
b = zeros(n,1); 
for i = 1:n 
    x_i = i*h; 
    if ((norm([x_i] - [0.3])<= 0.1) || (norm([x_i] - [0.6]) <= 0.1))
        b(i) = 1; 
    end
end
end 
%2D laplacian
function [A,b] = Lap2D(n) 
h = 1/(n+1); 
e = ones(n,1); 
% We form the 1D stencil matrices in each dimension. 
Ax = spdiags([e -2*e e], -1:1, n,n); 
Ay = spdiags([e -2*e e], -1:1, n,n); 
A2D = kron(Ax,Ay); 
% Build the left-hand-side matrix of finite difference: 
A = -(1/(h^2))*A2D; 
% Store as sparse matrix: 
A = sparse(A); 
% Build right-hand-side matrix. 
b = zeros(n,n); 
for i = 1:n 
    for j = 1:n 
        x_i = i*h; 
        y_j = j*h; 
        if ((norm([x_i,y_j] - [0.3,0.3])<= 0.1) || (norm([x_i, y_j] - [0.6,0.6]) <= 0.1))
            b(i,j) = 1; 
        end
    end 
end
% flatten to column matrix: 
b = reshape(b, [size(b,1)*size(b,2), 1]); 
end 
% 3D laplacian
function [A,b] = Lap3D(n)

h = 1/(n+1); 
e = ones(n,1); % ones element vector. 
% the sparse matrix in each dimension x,y,z: 
Ax = spdiags([e -2*e e], -1:1, n,n);
Ay = Ax; 
Az = Ax; 
% the identity matrix in each dimension x,y,z: 
Ix = speye(n); 
Iy = Ix; 
Iz = Ix; 
% build 3D matrix: 
A3D = kron(kron(Iz,Iy),Ax) + kron(kron(Iz,Ay), Ix) + kron(kron(Az,Iy), Ix);
% Build the LHS of finite diffrence: 
A = -(1/(h^2))*A3D; 
% Store as sparse matrix: 
A = sparse(A); 
% Build RHS matrix: 
b = zeros(n,n,n); 
for i = 1:n 
    for j = 1:n 
        for k = 1:n 
            x_i = i*h; 
            y_j = j*h; 
            z_k = k*h; 
            if ((norm([x_i,y_j,z_k] - [0.3,0.3,0.3])<= 0.1) || (norm([x_i, y_j,z_k] - [0.6,0.6,0.6]) <= 0.1))
                b(i,j,k) = 1; 
            end 
        end 
    end 
end 
% flatten to 1D matrix: 
b = reshape(b, [size(b,1)*size(b,2)*size(b,3),1]); 
end 