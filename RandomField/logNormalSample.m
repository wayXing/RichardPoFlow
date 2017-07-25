function y=logNormalSample(mu_Y,Sigma_Y)
% Samples from a log-normal distribution log(Y)=X, where X~N(mu_X,Sigma_X).
%
%
% Input parameters:
%   mu_Y            -[d * 1] vector of mean of Y
%   Sigma_Y         -[d * d] matrix of covariance of Y
%
% Output parameters:
%   y               -[n * 1] samples
%
% Examples: 
%
% See also: 
%
% Author:   Wei Xing
% History:  24/07/2017  file created
%
%% Main
%
% calculate corresponding mu_X, Sigma_X.
% Todo: speed could be improved due to symmetricity of Sigma 
Sigma_X=log(Sigma_Y./(mu_Y*mu_Y')+ 1);
mu_X=log(mu_Y)-diag(Sigma_X)./2;

%sampling form X
% y = exp( mvnrnd( mu_X , Sigma_X , Simulations ));
y = exp( mvnrnd( mu_X , Sigma_X  ));


