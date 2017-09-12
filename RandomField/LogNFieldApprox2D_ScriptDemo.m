% function []=LogNFieldApprox2D_ScriptDemo()
% Log-normal random field script (non functions calling) demonstration in
% 2D
% Approximation including: Kl approximation 
%                          interpolation for high resolution
% Equation:
% Y ~logNormal(mu_Y,sigma_Y)
% sigma_Y_ij=exp(-|x_i-x_j|/c)*sigma
%
% Log: 
%
% Author:   Wei Xing
% History:  10/09/2017  file created

clear

%% set parameter
%overall parameters 
lengthX=20;
deltaX=1;
nX=lengthX/deltaX+1;

lengthY=20;
deltaY=1;
nY=lengthY/deltaY+1;

seed=102;
nSample=1;
muY=10;
lengthScale=5;
Deviation=0.1*muY;

% Kl parameters
nKl=50;

% interpolation parameters
fineDeltaX=0.1;
fineNX=lengthX/fineDeltaX+1;
fineDeltaY=0.1;
fineNY=lengthY/fineDeltaY+1;
fineD=fineNX*fineNY;



%% Main
[X,Y] = ndgrid(0:deltaX:lengthX,0:deltaY:lengthY);
location=[X(:),Y(:)];

distance = pdist(location);
distanceMatrix = squareform(distance);

% Calculate covariance matrix of Y
% MODIFY for richer structure
coarseSigmaY=exp(-distanceMatrix./lengthScale) .*Deviation^2;    

[muX,SigmaX]=LogN2N(muY*ones(size(location,1),1),coarseSigmaY);

%% KL on X
% nKl=50;
[klBasisX,klEigenValueX,~] = svds(SigmaX,nKl); 

energyX=diag(klEigenValueX);
energyRatioX = cumsum(energyX)./sum(energyX);


sample= randn(nKl,nSample);
% x=klBasisX*sqrt(klEigenValueX)*sample+repmat(muX,1,nSample);

kTemp=klBasisX*sqrt(klEigenValueX)*diag(sample);
kTemp=cumsum(kTemp,2)+repmat(muX,1,nKl);

K = reshape(exp(kTemp),nX,nY,nKl);

%% Interpolation on basis of coarse grid

[coarseX,coarseY] = ndgrid(0:deltaX*2:lengthX,0:deltaY*2:lengthY); %have to be just fine times

coarseLocation=[coarseX(:),coarseY(:)];

coarseDistance = pdist(coarseLocation);
CoarseDistanceMatrix = squareform(coarseDistance);

% Calculate covariance matrix of Y
% MODIFY for richer structure
coarseSigmaY=exp(-CoarseDistanceMatrix./lengthScale) .*Deviation^2;    

[coarseMuX,coarseSigmaX]=LogN2N(muY*ones(size(coarseLocation,1),1),coarseSigmaY);
[coarseKlBasis,coarseklEigenValue,~] = svds(coarseSigmaX,nKl); 

for i=1:nKl
%     surface=scatteredInterpolant(coarseLocation,klBasis(:,i),'natural');
%     surface=scatteredInterpolant(coarseLocation,klBasis(:,i),'nearest');
    surface=scatteredInterpolant(coarseLocation,coarseKlBasis(:,i));
    interpolatedKlBasis(:,i)=surface(location);
end


% interpolatedMuX=log(muY)-diag(fineKlBasis*klEigenValue*fineKlBasis')./2; %not feasible for high resolution
interpolatedMuX=repmat(muX(1,1),size(interpolatedKlBasis,1),1);    %this is wrong solution
% 
%         kTemp=interpolatedKlBasis*sqrt(klEigenValue)*sample+repmat(interpolatedMuX,1,nSample);
kTemp=interpolatedKlBasis*sqrt(coarseklEigenValue)*diag(sample);
kTemp=cumsum(kTemp,2)+repmat(interpolatedMuX,1,nKl);

interpolatedK = reshape(exp(kTemp),nX,nY,nKl);
        


%% Plot
%     figure(1)
%     set(gca, 'YScale', 'log')
%     plot(energyRatioX,'-r')
%     hold on
%     plot(energyRatioY,'.-k')
%     hold off
%     legend('Energy of X', 'Energy of Y')
%     title(sprintf('Energy Ratio'))
%     % set(gca, 'YScale', 'log')

showIndex=1:2;
figure(1)
basis=klBasisX*klEigenValueX;
plot(basis(:,showIndex),'-')
hold on 
basis=interpolatedKlBasis*coarseklEigenValue;
plot(basis(:,showIndex),'--')
hold off 
title(sprintf('Eigen function compare'))



figure(2)
pcolor(X,Y,K(:,:,end))
shading interp;
colormap jet;
colorbar
title(sprintf('True field'))


for i=1:1:nKl
    figure(3)
    subplot(1,2,1)
    pcolor(X,Y,K(:,:,i))
    shading interp;
    colormap jet;
    colorbar
    title(sprintf('nKL=%d, energy=%f',i,energyRatioX(i)))
    
    subplot(1,2,2)
    pcolor(X,Y,interpolatedK(:,:,i))
    shading interp;
    colormap jet;
    colorbar
    title(sprintf('upscaleK nKL=%d, energy=%f',i,energyRatioX(i)))
    
    
    
%     contourf(X,Y,K(:,:,i))
% %     colormap(hot)
%     shading interp;
%     colorbar
%     title(sprintf('Permeability field'))
    frame(i)=getframe;
end