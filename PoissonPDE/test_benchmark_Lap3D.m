function test_benchmark_Lap3D
clc; 
clear all;
close all; 

m = [20,30,40,50]; 
N_size = m.^3; % problem size of matrix A. 
dim = '3D';
maxit = 2000;
tol = 1e-12;

for i=1:length(m)
   m_grids = m(i); 
   [A,b] = Laplacian(m_grids, dim); 
   problem_size = size(A,1); 
   u0 = zeros(problem_size,1); % generate initial guess of solution.
   fprintf("\r\n---Problem Size N=%d---\r\n", problem_size);
   
   
   fprintf("\r\n----GMRES----\r\n"); 
   method = "GMRES"; 
   [~,~,relres_gmres, iter_gmres, resvec_gmres,elapsed_time] = KrylovMethod(A,b,u0,maxit,tol,method); 
   fprintf("\r\nThe GMRES method converges in %d iterations with relative residual %.2g in %.3f (s)\r\n", iter_gmres, resvec_gmres(end)/norm(b),elapsed_time);
   itersGMRES(i) = iter_gmres; 
   resval_gmres(i) = resvec_gmres(end); 
   
   
   
   fprintf("\r\n----BiCG----\r\n"); 
   method = "BiCG"; 
   [~,~,relres_bicg, iter_bicg,resvec_bicg, elapsed_time] = KrylovMethod(A,b,u0,maxit,tol,method); 
   fprintf("\r\nThe BiCG method converges in %d iterations with relative residual %.2g in %.3f (s) \r\n", iter_bicg, resvec_bicg(end)/norm(b),elapsed_time);
   itersBicg(i) = iter_bicg;
   resval_bicg(i) = resvec_bicg(end); 
   
   fprintf("\r\n----IDR(s)----\r\n"); 
   method = "IDRs"; 
   
   s4 = 4;
   [~,~,relres_idr4, iter_idr4,resvec_idr4, elapsed_time] = KrylovMethod(A,b,u0,maxit,tol,method,s4); 
   fprintf("\r\nThe IDR(s) with s=%d method converges in %d iterations with relative residual %.2g in %.3f (s)\r\n",s4,iter_idr4, resvec_idr4(end)/norm(b),elapsed_time);
   itersIDR4(i) = iter_idr4; 
   resval_idr4(i) = resvec_idr4(end); 
   
    s8 = 8;
   [~,~,relres_idr8, iter_idr8,resvec_idr8,elapsed_time] = KrylovMethod(A,b,u0,maxit,tol,method,s8); 
   fprintf("\r\nThe IDR(s) with s=%d method converges in %d iterations with relative residual %.2g in %.3f (s) \r\n",s8,iter_idr8, resvec_idr8(end)/norm(b),elapsed_time);
   itersIDR8(i) = iter_idr8; 
   resval_idr8(i) = resvec_idr8(end); 
   
   
    s100 = 100;
   [~,~,relres_idr100, iter_idr100,resvec_idr100,elapsed_time] = KrylovMethod(A,b,u0,maxit,tol,method,s100); 
   fprintf("\r\nThe IDR(s) with s=%d method converges in %d iterations with relative residual %.2g in %.3f (s) \r\n",s100,iter_idr100, resvec_idr100(end)/norm(b),elapsed_time);
   itersIDR100(i) = iter_idr100; 
   resval_idr100(i) = resvec_idr100(end); 
   
   
   if (m_grids == m(end))
      f1 = figure(1); 
      f1.Position = [100 100 1200 1000];
      
      resvec_gmres = resvec_gmres./norm(b); 
      resvec_bicg = resvec_bicg./norm(b); 
      resvec_idr4 = resvec_idr4./norm(b); 
      resvec_idr8 = resvec_idr8./norm(b); 
      resvec_idr100 = resvec_idr100./norm(b); 
      
      iters_gmres = 1:iter_gmres; 
      iters_bicg = 1:iter_bicg; 
      iters_idr4 = 1:iter_idr4; 
      iters_idr8 = 1:iter_idr8;
      iters_idr100 = 1:iter_idr100;
      
      % GMRES: 
      semilogy(iters_gmres, resvec_gmres, '--m+', 'LineWidth',4, 'MarkerSize',12);
      hold on; 
      % BiCG:
      semilogy(iters_bicg, resvec_bicg, '--r<', 'LineWidth',4, 'MarkerSize',12); 
      % IDR(s): 
      semilogy(iters_idr4, resvec_idr4,'-.kd', 'LineWidth',4, 'MarkerSize',12,'MarkerFaceColor', 'black');
      semilogy(iters_idr8, resvec_idr8,'--cp', 'LineWidth',4, 'MarkerSize',12,'MarkerFaceColor', 'cyan');
      semilogy(iters_idr100, resvec_idr100,'-.h', 'LineWidth',4, 'MarkerSize',12, 'Color', '#D95319', 'MarkerFaceColor', '#D95319');
      hold off;  
      grid on; 
      grid minor; 
      legend('gmres','bicg','idr(4)','idr(8)','idr(100)','Location','best', 'FontSize',20);
      xlabel('k iteration(s)', 'Interpreter', 'Latex','FontSize',35); 
      ylabel('$\log(\frac{||r_{k}||}{||b||} )$ ', 'Interpreter', 'Latex', 'FontSize',35);
      titleText = sprintf('3D Poisson with $N=%d$', problem_size); 
      title(titleText,'Interpreter', 'Latex', 'FontSize',40); 
      ax = gca; 
      ax.XAxis.FontSize = 30; 
      ax.YAxis.FontSize = 30;
      
   end    
end


f2 = figure(2); 
f2.Position = [100 100 800 500]; 
% GMRES: 
plot(N_size,itersGMRES,'-mo', 'LineWidth', 5, 'MarkerSize', 20, 'MarkerFaceColor', 'm');
hold on; 
plot(N_size, itersBicg,'-r<', 'LineWidth', 5, 'MarkerSize', 20, 'MarkerFaceColor', 'r'); 
plot(N_size, itersIDR4,'-.kd', 'LineWidth', 5, 'MarkerSize', 20,'MarkerFaceColor', 'black');
plot(N_size, itersIDR8,'-.cp', 'LineWidth', 5, 'MarkerSize', 20,'MarkerFaceColor', 'cyan');
plot(N_size, itersIDR100,'-.h', 'LineWidth',5, 'MarkerSize', 20, 'Color', '#D95319', 'MarkerFaceColor', '#D95319');
hold off; 
grid on; 
legend('gmres','bicg','idr(4)','idr(8)','idr(100)','Location', 'best','FontSize', 20);
xlabel('$log(N)$', 'Interpreter','Latex', 'FontSize', 30); 
ylabel('$log(k)$', 'Interpreter','Latex','FontSize', 30); 
titleText_f2 = sprintf('Loglog of Krylov methods in 3D problem');
title(titleText_f2, 'Interpreter', 'Latex','FontSize', 30); 
gca_f2 = get(f2, 'CurrentAxes'); 
set(gca_f2, 'YScale','log', 'XScale', 'log'); 
gca_f2.XAxis.FontSize = 25; 
gca_f2.YAxis.FontSize = 25;

end 