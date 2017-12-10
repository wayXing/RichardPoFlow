function y=permeaFieldApprox(location,lengthScale,muY,DeviationRatio,nSample,nKL)
% permeability field approximation generator (log-normal) given coordinates by Grid and measure lengthscale.
%
% Equation:
% Ks=Y ~logNormal(mu_Y,sigma_Y)
% sigma_Y_ij=exp(-|x_i-x_j|/c)
%
% Input parameters:
%   location        -[n * d] matrix with n data points of d-dimensional spatial location each row
%   lengthcale      -[1 * 1] scalar Correlation length. See equation. larger number means less stochastic field. Thus less smooth.
%   muY             -[1 * 1] scalar mean of Y 
%   DeviationRatio  -[1 * 1] scalar Deviation ratio
%   nSample         -[1 * 1] scalar number of sample
%   nKL             -[1 * 1] integer scalar <=n;
% Output parameters:
%   Ks              -[n * 1] matrix
%
% Examples: see Demo
%
% % Short description of example, followed by Matlab code line
% >> matlab code of example
%
% See also: 
% Author:   Wei Xing
% History:  24/07/2017  file created
%
%% Initial 
seed=101;
%
%% Main
%calculate distance matrix
distance = pdist(location);
distanceMatrix = squareform(distance);

% Calculate covariance matrix of Y
% TODO: MODIFY for richer structure
SigmaY=exp(-distanceMatrix./lengthScale) .*(muY*DeviationRatio)^2;      


% Conver to X covariance matrix and mean
SigmaX=log(SigmaY./(muY*muY')+ 1);
muX=log(muY)-diag(SigmaX)./2;


% KL decomposition of X; Normalization of X;
[klBasis,klEigenValue,~] = svds(SigmaX,nKL); 

%Generate independent normal samples 
% rng(seed);  %pseudo random

sample= randn(nKL,nSample);
x=klBasis*sqrt(klEigenValue)*sample+repmat(muX,1,nSample);

% a multi-variate log (multi) normal permeability field
y=exp(x);



end