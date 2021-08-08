% implement of denoise script
% To record the total iteration and cpu times for each of the linear 
% iterative method solver.
function Denoise
clear all; 
clc; 
close all; 
warning('off');
K = 10; 
tol = 1e-3; 
maxiter = 10000; 
m_grid = [128]; %[64];  
alpha = 8e-3;
iterLoop = size(m_grid,2);
windowSize = [2000, 1500]; % the size of figure window. 
for i=1:iterLoop
    m = m_grid(i); 
    fprintf("Grid resolution %d x %d:\n", m,m); 
    fprintf("The alpha value: %e\r\n", alpha); 
    [u_exact, u_noisy]=set_image(m); 
    u0=FormRHS(u_noisy); % the noisy image (RHS)
    method_1 = 'GMRES';
    method_2 = 'BiCG';
    method_3 = 'IDRs'; 
    s1 = 4;
    s2 = 8; 
    s3 = 100;
    
    % GMRES method: 
    fprintf("GMRES: running \r\n"); 
    [u_gmres,iter_gmres,time_gmres]=DenoiseFunc(u0,K,alpha,maxiter,tol,method_1); 
    fprintf("GMRES method is finised in %d iters, %.3f (secs) \r\n", iter_gmres, time_gmres); 
    
    % BiCG method. 
    fprintf("BiCG: running \r\n"); 
    [u_bicg,iter_bicg,time_bicg]=DenoiseFunc(u0,K,alpha,maxiter,tol,method_2); 
    fprintf("BiCG method is finised in %d iters, %.3f (secs) \r\n", iter_bicg, time_bicg);
    
    % IDR(4) method. 
    fprintf("IDR(4): running\r\n"); 
    [u_idr4, iter_idr4, time_idr4] = DenoiseFunc(u0,K,alpha,maxiter,tol,method_3,s1);
    fprintf("IDR(4) method is finished in %d iters, %.3f (secs) \r\n", iter_idr4, time_idr4); 
    
    % IDR(8) method. 
    fprintf("IDR(8): running\r\n"); 
    [u_idr8, iter_idr8, time_idr8] = DenoiseFunc(u0,K,alpha,maxiter,tol,method_3,s2);
    fprintf("IDR(4) method is finished in %d iters, %.3f (secs) \r\n", iter_idr4, time_idr4);
    
    % IDR(100) method. 
    fprintf("IDR(100): running\r\n"); 
    [u_idr100, iter_idr100, time_idr100] = DenoiseFunc(u0,K,alpha,maxiter,tol,method_3,s3);
    fprintf("IDR(100) method is finished in %d iters, %.3f (secs) \r\n", iter_idr4, time_idr4);
    
    
    % reshape the denoised image into 2D image and plot with imagesc gray
    % scale.
    set(gcf,'Position',[200 200 windowSize(1) windowSize(2)]);
    figure(1);
    
    
    subplot(2,3,1); 
    imagesc(u_exact); 
    colormap(gray); 
    title('Exact Image');
    
    subplot(2,3,2); 
    u_denoised_gmres = reshape(u_gmres, [m m])'; 
    imagesc(u_denoised_gmres); 
    colormap(gray); 
    title_text = sprintf('GMRES: %d iters, %.3f (secs)', iter_gmres,time_gmres); 
    title(title_text); 
    
    
    subplot(2,3,3); 
    u_denoised_bicg = reshape(u_bicg, [m m])'; 
    imagesc(u_denoised_bicg); 
    colormap(gray); 
    title_text = sprintf('BiCG: %d iters, %.3f (secs)', iter_bicg,time_bicg); 
    title(title_text); 
    
    
    subplot(2,3,4); 
    u_denoised_idr4 = reshape(u_idr4, [m m])'; 
    imagesc(u_denoised_idr4); 
    colormap(gray); 
    title_text = sprintf('IDR(4): %d iters, %.3f (secs)', iter_idr4,time_idr4); 
    title(title_text); 
    
    
    
    subplot(2,3,5); 
    u_denoised_idr8 = reshape(u_idr8, [m m])'; 
    imagesc(u_denoised_idr8); 
    colormap(gray); 
    title_text = sprintf('IDR(8): %d iters, %.3f (secs)', iter_idr8,time_idr8); 
    title(title_text);
    
    
    subplot(2,3,6); 
    u_denoised_idr100 = reshape(u_idr100, [m m])'; 
    imagesc(u_denoised_idr100); 
    colormap(gray); 
    title_text = sprintf('IDR(100): %d iters, %.3f (secs)', iter_idr100,time_idr100); 
    title(title_text);
    
    f1 = 'denoised_img_';
    f2 = num2str(m); 
    filename = [f1 f2];
    saveas(gcf,filename,'png');  
end 

end