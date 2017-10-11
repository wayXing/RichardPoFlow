% function []=LogN2N_Demo
% Test LogN2N and N2LogN function with a 2D random field whose covariance
% is isometric

clear

%% set parameter
lengthX=40;
deltaX=4;
nX=lengthX/deltaX+1;

lengthY=40;
deltaY=4;
nY=lengthY/deltaY+1;

muY=10;
lengthcale=10;
DeviationRatio=0.4;


%% Main
d=nX*nY;
[X,Y] = ndgrid(0:deltaX:lengthX,0:deltaY:lengthY);
location=[X(:),Y(:)];



distance = pdist(location);
distanceMatrix = squareform(distance);

% Calculate covariance matrix of Y
% MODIFY for richer structure
SigmaY=exp(-distanceMatrix./lengthcale) .*(muY*DeviationRatio)^2;    

[muX,SigmaX]=LogN2N(muY*ones(d,1),SigmaY);

[muY2,SigmaY2]=N2LogN(muX,SigmaX);

%show different of original and map-remap result 
ErrorMu=sum(muY-muY2)
ErrorSigma=sum(SigmaY(:)-SigmaY2(:))


