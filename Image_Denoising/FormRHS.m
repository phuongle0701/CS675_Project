% we form the right-hand-side matrix in eq.(2) 
% Input: noisy image X. 
% Output: vector representation of noisy image X. 
function [u0] = FormRHS(X)
X_t=X'; % transpose X. 
u0 = X_t(:); % vectorize. 
end 