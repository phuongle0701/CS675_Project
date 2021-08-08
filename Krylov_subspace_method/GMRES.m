% function to implement GMRES algorithm
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
function [x,flag,relres,iter,resvec] = GMRES(A,b,x0,tol,maxit)
% Generalized Minimum Residual Method with restarting 
% Correspond to Algorithm 6.9  in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)"
r = b - A*x0; 
residual  = norm(r,2); 
iter = 1;
if (residual < tol)
    flag = 0; 
    relres = residual/norm(b,2); 
    resvec(iter) = residual; 
end 

n = size(A, 1);
restart  = min(n, 10); 
x = zeros(n, 1);
r = b - A * x;

residual = norm(r);
iter = 1;
resvec = zeros(maxit * restart + 1, 1);

relres = residual/norm(b,2);
resvec(iter) = residual;
flag = 1; 
% main loop: 
while ((iter < maxit) && (residual > tol))
    
	[V, H, beta] = Arnoldi_MGS(A, r, restart);
	[y, resnorms]  = UpperHessenLeastSquare(H, beta); 
	z = V(:, 1 : restart) * y(1 : restart);

	x = x + z;
	r = b - A * x;
    residual = norm(r,2);
	for j = 1 : restart
		iter = iter + 1;
		resvec(iter) = resnorms(j);
		residual = min(residual, resnorms(j));
		if (residual < tol)
            flag = 0;
            relres = residual/norm(b,2);
			break;
		end
	end
end

% column vector; 
resvec = resvec(1:iter);
end


