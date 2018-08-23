function [AIC] = aic(Residuals,n,k)

% Estimation of the Akaike Information Criteria (AIC)
% Lower values theoretically indicate better model
%
% Information criteria provide a tradeoff between goodness-of-fit 
% (which lowers the sum of squared residuals) and the model's complexity, 
% which is measured by the number of parameters k+1.
% sum of squared residuals applicable for least-squares fits
% 
% Information criteria have to be minimized over choices of k+1.
%
% Residuals - differences between fit and data 
% n - number of data points
% k - number of parameters (excluding variance)
%
% pass in results from multiple models
% 
% Since the penalty for additional regressors is larger in BIC, 
% this criterion tends to favor more parsimonious models than AIC.
%
% Version 1.0
% ---------------------------------------------------------------------------
%
% Copyright Notice:
%
% You are allowed to use the code for your private and commercial purpose
% and change it to meet your requirements. You are not allowed to
% redistribute or sell it as a whole or fragments of it. When using it,
% cite it.
% 
% Copyright 2011 | Lon Bueckins | leon.bueckins@googlemail.com
%
% If you have any questions or suggestions for improvements, feel free to
% contact me.
%
%

% ideally n is large compared to k...

AIC = log(1./n * (Residuals'*Residuals)) + (2.*(k+1))./n;  % normalized, with k+1 to account for estimating variance

AIC = n .* log((Residuals'*Residuals) ./ n) + 2*(k+1); % raw

AIC = n .* log((Residuals'*Residuals)./n) + 2*(k+1) .* (n ./ (n-k-2) );  % small sample size corrected

% Diff = AIC - min(AIC);
% W = exp(-Diff/2) / sum(exp(-Diff/2));  
end