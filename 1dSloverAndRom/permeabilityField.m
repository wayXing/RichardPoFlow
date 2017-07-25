function Ks=permeabilityField(X,lengthcale)
% permeability generator (log-normal) given X coordinates and measure lengthscale.
%
% Equation:
% log(ks)=z~N(0,cov(x1,x2))
% cov[z_{12}]=cov(x1,x2)=exp(-|x1-x2|/c) and E[z]=0
%
% Input parameters:
%   X               -[n * d] matrix with n data points for spatial location each row
%   lengthcale      -larger number means less stochastic field. Thus less smooth.
% Output parameters:
%   Ks              -[n * 1] matrix
%   type
%
% Examples: see Demo
%
% % Short description of example, followed by Matlab code line
% >> matlab code of example
%
% See also: 
% Author:   Wei Xing
% History:  10/05/2017  file created


[nX,dimX]=size(X);

%calculate distance matrix
distance = pdist(X);
distanceMatrix = squareform(distance);

%calculate covariance matrix
covMatrix=exp(-distanceMatrix./lengthcale);    

% KL decomposition on covariance matrix via SVD/eigen decomposition
% [klBasis,klEigenValue] = eigs(covMatrix,nY*nX); 
[klBasis,klEigenValue,~] = svds(covMatrix,nX); 


% [nKlBasis,~]=sizes(klBasis);


%Generate independent normal samples 
seed=101;
rng(seed);
sample= randn(nX,1);

%make multivariate Gaussian distributions with samples. zero mean.
%Covariance specified though KL basis.
Ks=klBasis*sqrt(klEigenValue)*sample;
% Ks=reshape(Ks,nY,nX);

% a log (multi) normal permeability field
Ks=exp(Ks);

end