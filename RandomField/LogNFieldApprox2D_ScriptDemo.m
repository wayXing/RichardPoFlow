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
%           09/10/2017  Modify. Add upsale and test on it 

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
nKl=10;

% interpolation parameters
fineDeltaX=0.5;
fineNX=lengthX/fineDeltaX+1;
fineDeltaY=0.5;
fineNY=lengthY/fineDeltaY+1;
fineD=fineNX*fineNY;



%% Coarse grid 
[coarseX,coarseY] = ndgrid(0:deltaX:lengthX,0:deltaY:lengthY);
% [coarseX,coarseY] = meshgrid(0:deltaX:lengthX,0:deltaY:lengthY);

coarseLocations=[coarseX(:),coarseY(:)];

coarseDistances = pdist(coarseLocations);
coarseDistanceMatrix = squareform(coarseDistances);

% Calculate covariance matrix of Y
% MODIFY for richer structure
coarseSigmaY=exp(-coarseDistanceMatrix./lengthScale) .*Deviation^2;    

[coarseMuX,coarseSigmaX]=LogN2N(muY*ones(size(coarseLocations,1),1),coarseSigmaY);


%% Fine grid 
[fineX,fineY] = ndgrid(0:fineDeltaX:lengthX,0:fineDeltaY:lengthY);
fineLocations=[fineX(:),fineY(:)];

fineDistances = pdist(fineLocations);
fineDistanceMatrix = squareform(fineDistances);

% Calculate covariance matrix of Y
% MODIFY for richer structure
fineSigmaY=exp(-fineDistanceMatrix./lengthScale) .*Deviation^2;    

[fineMuX,fineSigmaX]=LogN2N(muY*ones(size(fineLocations,1),1),fineSigmaY);

%% KL on fine grid
% nKl=50;
[klBasisFineX,klEigenValueFineX,~] = svds(fineSigmaX,nKl); 

energyX=diag(klEigenValueX);
energyRatioX = cumsum(energyX)./sum(energyX);


sample= randn(nKl,nSample);
% x=klBasisX*sqrt(klEigenValueX)*sample+repmat(muX,1,nSample);

kTemp=klBasisX*sqrt(klEigenValueX)*diag(sample);
kTemp=cumsum(kTemp,2)+repmat(muX,1,nKl);

K = reshape(exp(kTemp),nX,nY,nKl);


%% KL on coarse grid














%% Interpolation on basis of coarse grid






% [coarseX,coarseY] = ndgrid(0:deltaX:lengthX,0:deltaY:lengthY); %have to be just fine times

% coarseLocation=[coarseX(:),coarseY(:)];

% coarseDistances = pdist(coarseLocation);
% CoarseDistanceMatrix = squareform(coarseDistances);

% Calculate covariance matrix of Y
% MODIFY for richer structure
% coarseSigmaY=exp(-CoarseDistanceMatrix./lengthScale) .*Deviation^2;    

% [coarseMuX,coarseSigmaX]=LogN2N(muY*ones(size(coarseLocation,1),1),coarseSigmaY);
[coarseKlBasis,coarseklEigenValue,~] = svds(coarseSigmaX,nKl); 

iMethod=2;
switch iMethod 
    case 1
        for i=1:nKl
        %     surface=scatteredInterpolant(coarseLocation,klBasis(:,i),'natural');
        %     surface=scatteredInterpolant(coarseLocation,klBasis(:,i),'nearest');
            surface=scatteredInterpolant(coarseLocation,coarseKlBasis(:,i));
            interpolatedKlBasis(:,i)=surface(coarseLocations);
        end

    case 2
        for i=1:nKl
            iCoarseKlBasisField = reshape(coarseKlBasis(:,i),nX,nY);
%             iInterpolatedKlBasisField = interp2(coarseLocation(:,1),coarseLocation(:,2),iCoarseKlBasisField,...
%                                                 location(:,1),location(:,2));
                                            
            iInterpolatedKlBasisField = interp2(coarseX,coarseY,iCoarseKlBasisField,...
                                                fineX,fineY,'spline');
                                            
            interpolatedKlBasis(:,i)= iInterpolatedKlBasisField(:);                               
        end
        
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
pcolor(coarseX,coarseY,K(:,:,end))
shading interp;
colormap jet;
colorbar
title(sprintf('True field'))


for i=1:1:nKl
    figure(3)
    subplot(1,2,1)
    pcolor(coarseX,coarseY,K(:,:,i))
    shading interp;
    colormap jet;
    colorbar
    title(sprintf('nKL=%d, energy=%f',i,energyRatioX(i)))
    
    subplot(1,2,2)
    pcolor(coarseX,coarseY,interpolatedK(:,:,i))
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