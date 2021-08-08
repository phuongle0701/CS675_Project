% script to combine the three methods to solve the system: Ax = b
% Output: 
% x - approximate solution 
% flag - converged or not. 
% relres - relative residual 
% iter  - converged iteration
% resvec - vector storing the residual norms || b - Ax || 
% elapsed_time: time it take to run through the iterative methods. 
function [x,flag,relres,iter, resvec, elapsed_time] = KrylovMethod(A,b,x0,maxit,tol,method,s)

if (nargin > 6)
    sVal = s; % use IDR(s)
else
    sVal = 0; % not use IDR(s). 
end

switch method
    case 'GMRES'
        tic;
        [x,flag,relres,iter,resvec] = gmres(A,b,[],tol,maxit,[],[],x0); 
        elapsed_time = toc; 
        iter = iter(2) + 1; % iteration at 0.
    case 'BiCG'
        tic;
        [x,flag,relres,iter,resvec] = BiCG(A,b,x0,tol,maxit); 
        elapsed_time = toc; 
    case 'IDRs'
        tic; 
        [x,flag,relres,iter,resvec] = IDRs(A,b,sVal,tol,maxit,[],[],x0,0);
        elapsed_time = toc; 
        iter = iter + 1; % G_i space at 0. 
    % solve with built-in     
    case 'Backslash'
        tic; 
        x = A\b; 
        elapsed_time = toc; 
        flag = 0; 
        relres = []; 
        iter = []; 
        resvec = []; 
end 

end 