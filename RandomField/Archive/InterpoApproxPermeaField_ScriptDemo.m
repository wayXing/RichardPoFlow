% function y=InterpoApproxPermeaField_ScriptDemo()
% Interpolation + Approximate permeability generator (log-normal) using KL decomposition
%
% Equation:
% Y ~logNormal(mu_Y,sigma_Y)
% sigma_Y_ij=exp(-|x_i-x_j|/c)*sigma
%
% Log: -permeaField_ScriptDemo();       Generate random field and plot
%      -ApproxPermeaField_ScriptDemo(); Generate approximate random field and plot
%
% Author:   Wei Xing
% History:  10/09/2017  file created

clear

%% set parameter
lengthX=40;
deltaX=1;
nX=lengthX/deltaX+1;

lengthY=40;
deltaY=1;
nY=lengthY/deltaY+1;

fineDeltaX=0.1;
fineNX=lengthX/fineDeltaX+1;
fineDeltaY=0.1;
fineNY=lengthY/fineDeltaY+1;
fineD=fineNX*fineNY;

seed=102;
nSample=1;
muY=10;
lengthcale=5;
Deviation=0.01*muY;

%% Main
d=nX*nY;        %Dimension of random vector
[X,Y] = ndgrid(0:deltaX:lengthX,0:deltaY:lengthY);
location=[X(:),Y(:)];

distance = pdist(location);
distanceMatrix = squareform(distance);

% Calculate covariance matrix of Y
% MODIFY for richer structure
SigmaY=exp(-distanceMatrix./lengthcale) .*Deviation^2;    

[muX,SigmaX]=LogN2N(muY*ones(d,1),SigmaY);

%% KL on X
nKl=50;

[klBasisX,klEigenValueX,~] = svds(SigmaX,nKl); 

energyX=diag(klEigenValueX);
energyRatioX = cumsum(energyX)./sum(energyX);


sample= randn(nKl,nSample);
% x=klBasisX*sqrt(klEigenValueX)*sample+repmat(muX,1,nSample);

kTemp=klBasisX*sqrt(klEigenValueX)*diag(sample);
kTemp=cumsum(kTemp,2)+repmat(muX,1,nKl);

K = reshape(exp(kTemp),nX,nY,nKl);

%% Interpolation 
[fineX,fineY] = ndgrid(0:fineDeltaX:lengthX,0:fineDeltaY:lengthY);
location=[X(:),Y(:)];

% for i=1:size(klBasisX,2)
%     fineKlBasisX(:,i)=griddata(X(:),Y(:),klBasisX(:,i),fineX(:),fineY(:));
% end

for i=1:size(klBasisX,2)
    surface=scatteredInterpolant(X(:),Y(:),klBasisX(:,i));
    fineKlBasisX(:,i)=surface(fineX(:),fineY(:));
end

kTemp=fineKlBasisX*sqrt(klEigenValueX)*diag(sample);

kTemp=cumsum(kTemp,2)+repmat(muX(1,1),fineD,nKl);       %TODO generalize

interpoK = reshape(exp(kTemp),fineNX,fineNY,nKl);


%% KL on Y
%         [klBasisY,klEigenValueY,~] = svds(SigmaY,d); 
% 
%         energyY=diag(klEigenValueY);
%         energyRatioY = cumsum(energyY)./sum(energyY);



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

figure(2)
pcolor(X,Y,K(:,:,end))
shading interp;
colormap jet;
colorbar
title(sprintf('True field'))


for i=1:1:d
    figure(3)
    subplot(1,2,1)
    pcolor(X,Y,K(:,:,i))
    shading interp;
    colormap jet;
    colorbar
    title(sprintf('nKL=%d, energy=%f',i,energyRatioX(i)))
    
    subplot(1,2,2)
    pcolor(fineX,fineY,interpoK(:,:,i))
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