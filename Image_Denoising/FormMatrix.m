% To form the matrix corresponds to equation (2)
% Input: u (approximation solution vector); alpha (coefficient in eq.(2)). 
% Output: The matrix in eq.(2). 

function [A] = FormMatrix(u,alpha)
% sizes
n = size(u,1); 
m = sqrt(n); 
h = 1/(m+1); % grid size. 

% Hyper-parameter. 
beta = 1e-6; % 10^(-6)

% input reshape: 
u = reshape(u,[m,m])'; 

% Setup storage. 
AW = zeros(m,m); 
AE = zeros(m,m); 
AS = zeros(m,m); 
AN = zeros(m,m); 
AC = zeros(m,m); 

% helper function for stencil point increment/decrement. 
    function output=s(i,j)
       if (i<=m && i>= 1 && j<= m && j>=1)
          output = u(i,j); 
       else
          output = 0; 
       end 
    end 

% function to implement the equation inside each direction 
% output = 1 / (2 * sqrt(((a - b) / h) ^ 2 + ((c - d) / h) ^ 2 + beta));
% a,b,c,d: stencil points 
    function output = formula(a,b,c,d)
        output = 1/(2 * sqrt(((a-b)/h)^2 + ((c-d)/h)^2 + beta));
    end
% implement the generalize formula in each direction
% a,b,c,d,e,f,g,i: stencil points. 
    function output = formDirect(a,b,c,d,e,f,g,i)
        output = (-alpha/(h^2))*(formula(a,b,c,d) + formula(e,f,g,i));
    end 

% We form the direction matrix:
for i=1:m
    for j=1:m
        % West direction
      AW(i,j) = formDirect(s(i,j), s(i-1,j),s(i, j),s(i,j-1),s(i, j),s(i-1,j),s(i-1,j+1),s(i-1,j));
        % East direction
      AE(i,j) = formDirect(s(i+1,j),s(i,j),s(i+1,j),s(i+1,j-1),s(i+1,j),s(i,j),s(i,j+1),s(i,j));
        % South direction
      AS(i,j) = formDirect(s(i,j),s(i-1, j), s(i,j), s(i,j-1), s(i+1, j-1), s(i,j-1), s(i,j), s(i,j-1));
        % North direction
      AN(i,j) = formDirect(s(i+1,j),s(i,j),s(i,j+1),s(i,j),s(i,j+1),s(i-1,j+1),s(i,j+1),s(i,j));
        % Center
      AC(i,j) = - (AW(i,j) + AE(i,j) + AS(i,j) + AN(i,j))+1;
    end
end 

% vectorize: 
AC_t = AC';
AN_t = AN';
AS_t = AS';
AE_t = AE(1:m-1,:)';
AW_t = AW(2:m,:)';
AC_vec = AC_t(:);
AN_vec = AN_t(:);
AS_vec = AS_t(:);
AE_vec = AE_t(:);
AW_vec = AW_t(:);
AN_vec = AN_vec(1:n - 1);
AS_vec = AS_vec(2:n);

% matrix in equation (2). 
A = diag(AC_vec)+diag(AN_vec,1)+diag(AS_vec,-1)+diag(AE_vec,m)+diag(AW_vec,-m);

% assume zero boundary conditions  
for i=1:(m-1)
  A(i * m, i * m + 1) = 0;
  A(i * m + 1, i * m) = 0;
end

% Finalize the result matrix as sparse matrix for storage. 
A = sparse(A);

end
