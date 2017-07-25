% function y=ApproxPermeaField_ScriptDemo()
% Approximate permeability generator (log-normal) using KL decomposition
%
% Equation:
% Y ~logNormal(mu_Y,sigma_Y)
% sigma_Y_ij=exp(-|x_i-x_j|/c)*sigma
%
% Log: -permeaField_ScriptDemo();       Generate random field and plot
%      -ApproxPermeaField_ScriptDemo(); Generate approximate random field and plot
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
[klBasisX,klEigenValueX,~] = svds(SigmaX,d); 

energyX=diag(klEigenValueX);
energyRatioX = cumsum(energyX)./sum(energyX);


sample= randn(d,nSample);
% x=klBasisX*sqrt(klEigenValueX)*sample+repmat(muX,1,nSample);

kTemp=klBasisX*sqrt(klEigenValueX)*diag(sample);
kTemp=cumsum(kTemp,2)+repmat(muX,1,d);

K = reshape(exp(kTemp),nX,nY,d);

%% KL on Y
[klBasisY,klEigenValueY,~] = svds(SigmaY,d); 

energyY=diag(klEigenValueY);
energyRatioY = cumsum(energyY)./sum(energyY);



%% Plot
figure(1)
set(gca, 'YScale', 'log')
plot(energyRatioX,'-r')
hold on
plot(energyRatioY,'.-k')
hold off
legend('Energy of X', 'Energy of Y')
title(sprintf('Energy Ratio'))
% set(gca, 'YScale', 'log')

figure(2)
pcolor(X,Y,K(:,:,end))
shading interp;
colormap jet;
colorbar
title(sprintf('True field'))


for i=1:10:d
    figure(3)
    pcolor(X,Y,K(:,:,i))
    shading interp;
    colormap jet;
    colorbar
    title(sprintf('nKL=%d, energy=%f',i,energyRatioX(i)))
    
%     contourf(X,Y,K(:,:,i))
% %     colormap(hot)
%     shading interp;
%     colorbar
%     title(sprintf('Permeability field'))
    frame(i)=getframe;
end