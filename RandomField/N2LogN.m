function [muY,SigmaY]=N2LogN(muX,SigmaX)
% convert x normal distribution to exp(x) log normal distribution Y's mean
% and covariance 
%
% Equation:
%   X~N(muX,SigmaX) 
%   Y~LogNormal(muY,SigmaY) 
%   log(Y)=X
%    
% Input parameters:
%   muX             -[d * 1] mean of X 
%   SigmaX          -[d * d] Covariance matrix of X
% Output parameters:
%   muY             -[d * 1] mean of Y 
%   SigmaY          -[d * d] Covariance matrix of Y
%
% Examples: see Demo
% See also: 
%
% Author:   Wei Xing
% History:  25/07/2017  file created
%
%
%% Main
% Conver to Y covariance matrix and mean
muY   =exp(muX+diag(SigmaX)./2);
SigmaY= muY*muY' .* (exp(SigmaX)-1);



end