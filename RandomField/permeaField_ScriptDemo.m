% function y=permeaField_ScriptDemo()
% permeability generator (log-normal)
%
% Equation:
% Y ~logNormal(mu_Y,sigma_Y)
% sigma_Y_ij=exp(-|x_i-x_j|/c)
%
% Log: -permeaField_ScriptDemo() Generate random field and plot
%
% Author:   Wei Xing
% History:  25/07/2017  file created

clear

%% set parameter
lengthX=40;
deltaX=1;
nX=lengthX/deltaX+1;

lengthY=40;
deltaY=1;
nY=lengthY/deltaY+1;

seed=101;
nSample=10;
muY=10;
lengthcale=10;
DeviationRatio=0.05;
nKL=10;

%% Main
d=nX*nY;        %Dimension of random vector
[X,Y] = ndgrid(0:deltaX:lengthX,0:deltaY:lengthY);
location=[X(:),Y(:)];

distance = pdist(location);
distanceMatrix = squareform(distance);

% Calculate covariance matrix of Y
% MODIFY for richer structure
SigmaY=exp(-distanceMatrix./lengthcale) .*(muY*DeviationRatio)^2;    

[muX,SigmaX]=LogN2N(muY*ones(d,1),SigmaY);


% rng(seed);  %pseudo random
K = exp( mvnrnd( muX , SigmaX , nSample ))';
K = reshape(K,nX,nY,nSample);


for i=1:nSample
    figure(i)

    pcolor(X,Y,K(:,:,i))
    shading interp;
    colormap jet;
    colorbar
    title(sprintf('Permeability field'))
    
%     contourf(X,Y,K(:,:,i))
% %     colormap(hot)
%     shading interp;
%     colorbar
%     title(sprintf('Permeability field'))
    
end