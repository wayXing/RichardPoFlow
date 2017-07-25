function [muX,SigmaX]=LogN2N(muY,SigmaY)
% convert log normal distribution Y's mean and covariance to it
% correspondind generating x normal distribution.
%
% Equation:
%   X~N(muX,SigmaX) 
%   Y~LogNormal(muY,SigmaY) 
%   log(Y)=X
%
% Input parameters:
%   muY             -[d * 1] mean of Y 
%   SigmaY          -[d * d] Covariance matrix of Y
% Output parameters:
%   muX             -[d * 1] mean of X 
%   SigmaX          -[d * d] Covariance matrix of X
%
% Examples: see Demo
% See also: 
%
% Author:   Wei Xing
% History:  24/07/2017  file created
%
%
%% Main
% Conver to X covariance matrix and mean
SigmaX=log(SigmaY./(muY*muY')+ 1);
muX=log(muY)-diag(SigmaX)./2;



end