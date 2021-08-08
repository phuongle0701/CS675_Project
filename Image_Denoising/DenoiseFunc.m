% function to denoise the image with a Krylov space method
function [u_denoised, iter_sum,cpu_time] = DenoiseFunc(u0,K,alpha,maxit,tol,method,s)
%Setup
iter_sum = 0; 
u_image = u0;
% check the input parameters. 
if(nargin > 6)
   sVal = s; % use IDRs 
else 
    sVal = 4; % 
end 

% denoise the image: 
switch method
    case 'GMRES'
        tic; 
        for k=0:K
             A = FormMatrix(u_image,alpha);
             [u_image,~,~,iters,~,~] = KrylovMethod(A,u0,u_image,maxit,tol,method);
             iter_sum = iters + iter_sum;
        end 
        cpu_time = toc; 
        u_denoised = u_image; 
   
    case 'BiCG' 
        tic; 
        for k=0:K
            A = FormMatrix(u_image,alpha); 
            [u_image,~,~,iters,~,~] = KrylovMethod(A,u0,u_image,maxit,tol,method);
            iter_sum = iters+iter_sum;
        end
        cpu_time = toc; 
        u_denoised = u_image; 

    case 'IDRs'
        tic; 
        for k=0:K
            A = FormMatrix(u_image,alpha); 
            [u_image,~,~,iters,~,~] = KrylovMethod(A,u0,u_image,maxit,tol,method,sVal); 
            iter_sum = iters + iter_sum; 
        end 
        cpu_time = toc; 
        u_denoised = u_image; 

end 

end 





