% function to implement the BiCG algorithm in YSadd 
% Input: 
% Args: 
% A = matrix non-symmetric. 
% b = matrix of Ax = b
% x0 = initial guess of solution 
% tol = tolerance. 
% maxit = maximum number of iterations 

% Output: 
% x = approximate solution. 
% flag = return 0 if converged. 
% relres = final converged relative residuals (2-norm). 
% iter = number of iterations taken to converge. 
% resvec = vector to store the residuals at each iteration. 
function [x, flag,relres,iter,resvec] = BiCG(A,b,x0,tol, maxit)
n = size(A,1); 
x = x0; 
r = b - A*x; 
residual  = norm(r,2); 
iter = 1;
if (residual < tol)
    flag = 0; 
    relres = residual/norm(b,2); 
    resvec(iter) = residual; 
end 

resvec(1) = residual;

p = r; 
rs = r; 
ps = p; 
r2 = r' * rs;
AT = A'; 
Ap = A*p; 
ATps = AT * ps; 
flag = 1; 

while ((iter < maxit) && (residual >= tol))
    
    alpha = r2/(Ap' * ps); 
    x = x + alpha * p; 
    r = r - alpha * Ap; 
    rs = rs - alpha * ATps; 
    
    r2_old = r2; 
    r2 = r' * rs; 
    beta = r2/r2_old; 
    
    p = r + beta * p; 
    ps = rs + beta * ps; 
    Ap = A*p; 
    ATps = AT * ps; 
    
    
    iter = iter + 1; 
    residual = norm(r,2); 
    relres = residual/norm(b,2); 
    resvec(iter) = residual;   
end 

if (residual < tol) 
    flag = 0;
end 
% column vector: 
resvec = resvec(:);
end 