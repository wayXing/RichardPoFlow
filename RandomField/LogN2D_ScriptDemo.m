% function []=LogN2D_ScriptDemo()
% Log-normal random field script (non functions calling) demonstration in
% 2D using KL on fine and coarse grid. 
% Approximation including: Kl approximation 
%                          interpolation for high resolution
% Equation:
% Y ~logNormal(mu_Y,sigma_Y)
% sigma_Y_ij=exp(-|x_i-x_j|/c)*sigma
%
% Log: 
% WARMING: the code is not complete as the interpolated kl basis do not line up properly
%
% Author:   Wei Xing
% History:  09/10/2017  file created

clear; close all

%% set parameter
%TODO: Define covariance function here

seed=102;

nSample=5;
muY=10;
lengthScale=2;
sDeviation=0.1*muY;

% Kl parameters
nKl=50;

%% Grid1    a coarse grid 
% [grid1x,grid1y] = ndgrid(0:2:10,0:2:20);
[grid1x,grid1y] = meshgrid(0:1:10,0:1:10);


grid1Distances= pdist([grid1x(:),grid1y(:)]);
grid1DistanceMatrix = squareform(grid1Distances);

% Calculate covariance matrix of Y
SigmaY1=exp(-grid1DistanceMatrix./lengthScale) .*sDeviation^2;    
[MuX1,SigmaX1]=LogN2N(muY*ones(size(SigmaY1,1),1),SigmaY1);

%KL
[eigVec1,eigVal1,~] = svds(SigmaX1,nKl); 
klBasis1=eigVec1*sqrt(eigVal1);


%% Grid2    a fine grid 
% [grid2x,grid2y] = ndgrid(0:0.5:10,0:0.5:10);
[grid2x,grid2y] = meshgrid(0:0.5:10,0:0.5:10);


grid2Distances= pdist([grid2x(:),grid2y(:)]);
grid2DistanceMatrix = squareform(grid2Distances);

% Calculate covariance matrix of Y
SigmaY2=exp(-grid2DistanceMatrix./lengthScale) .*sDeviation^2;    
[MuX2,SigmaX2]=LogN2N(muY*ones(size(SigmaY2,1),1),SigmaY2);

%KL
[eigVec2,eigVal2,~] = svds(SigmaX2,nKl); 
klBasis2=eigVec2*sqrt(eigVal2);

%% TODO: line up kl basis of grid1&2 in the same direction 
% isVecDirectionSame=diag(sign(eigVec1(10,:).*klBasis1(10,:)),0);
% klBasis1= klBasis1* diag(sign(klBasis2(10,:).*klBasis1(10,:)),0);
% 
% % klBasis1= klBasis1* diag(sign(klBasis2*klBasis1',0);

%% Interpolation of fine grid (2) from coarse grid (1)

iMethod=3;
switch iMethod 
    case 1
        for i=1:nKl
        %     surface=scatteredInterpolant(coarseLocation,klBasis(:,i),'natural');
        %     surface=scatteredInterpolant(coarseLocation,klBasis(:,i),'nearest');
            surface=scatteredInterpolant(coarseLocation,grid(:,i));
            interpolatedKlBasis(:,i)=surface(coarseLocations);
        end

%     case 2  %TODO: interpolation the grid as a 1d vector 
%         for i=1:nKl
%             interpKlBasis2(:,i)= interp(,klBasis1(:,i),)
%            
%         end
        
        
    case 3  %interpolation as a 2D grid 
        
        MuXGrid1 = reshape(MuX1,size(grid1x,1),size(grid1x,2));                                   
        interpMuXGrid2 = interp2(grid1x,grid1y,MuXGrid1,...
                                                grid2x,grid2y,'cubic');
        
        for i=1:nKl
            klGrid1 = reshape(klBasis1(:,i),size(grid1x,1),size(grid1x,2));                                   
            interpKlGrid2 = interp2(grid1x,grid1y,klGrid1,...
                                                grid2x,grid2y,'cubic');
                                            
            interpKlBasis2(:,i)= interpKlGrid2(:);                               
        end
        
        % line up kl basis of klBasis2 and interpKlBasis2 in the same direction 
        isSameVecDirection=sign(diag(interpKlBasis2'*klBasis2));
        interpKlBasis2= interpKlBasis2* diag(isSameVecDirection); 
        
        
        
    case 4  %interpolation eigen vector as a 2D grid result
        for i=1:nKl
            iEigVecGrid1 = reshape(eigVec1(:,i),size(grid1x,1),size(grid1x,2));                                   
            iInterpEigVecGrid2 = interp2(grid1x,grid1y,iEigVecGrid1,...
                                                grid2x,grid2y,'spline');
                                            
            isVecDirectionSame=sign(iInterpEigVecGrid2(:)'*eigVec2(:,i));
                                
            interpKlBasis2(:,i)= iInterpEigVecGrid2(:)*sqrt(eigVal1(i,i))*isVecDirectionSame;
%                                  *sqrt(sqrt(size(grid2x,1)*size(grid2x,2))/sqrt(size(grid1x,1)*size(grid1x,2)));
%             interpKlBasis2(:,i)= iInterpEigVecGrid2(:)*sqrt(eigVal1(i,i))*isVecDirectionSame;
        end
        

end
interpSigmaX2=interpKlBasis2*interpKlBasis2';


%     %% line up kl basis of klBasis2 and interpKlBasis2 in the same direction 
%             isSameVecDirection=sign(diag(interpKlBasis2'*klBasis2));
%             interpKlBasis2= interpKlBasis2* diag(isSameVecDirection); 
% 
%     % klBasis1= klBasis1* diag(sign(klBasis2*klBasis1',0);

%% Re-order and re-direction  
BasisDist=pdist2(interpKlBasis2',klBasis2');
% BasisDist = squareform(basisDist);



%% Sampling
Theata= randn(nKl,nSample);

X2=klBasis2*Theata+repmat(MuX2,1,nSample);
Y2=exp(X2);

interpX2=interpKlBasis2*Theata+repmat(interpMuXGrid2(:),1,nSample);
interpY2=exp(interpX2);

%% visulizing covariance matrix
figure(1)
subplot(1,2,1)
imagesc(SigmaX2)
title(sprintf('Real Covariance'))

subplot(1,2,2)
imagesc(interpSigmaX2)
title(sprintf('interpolated Covariance'))


%% Plot kl basis 
showIndex=1:5;

figure(2)
plot(klBasis2(:,showIndex),'-')
hold on 
plot(interpKlBasis2(:,showIndex),'--')
hold off

%% Plot kl basis in grid
klGrid1 = reshape(klBasis1,size(grid1x,1),size(grid1x,2),[]); 
klGrid2 = reshape(klBasis2,size(grid2x,1),size(grid2x,2),[]); 
iInterpEigVecGrid2 = reshape(interpKlBasis2,size(grid2x,1),size(grid2x,2),[]); 

for i=1:1:10
    figure(10+i)
    
    subplot(1,3,1)
    surf(grid1x,grid1y,klGrid1(:,:,i))
    title(sprintf('Coarse KL %d-st KL',i))
    
    subplot(1,3,2)
    surf(grid2x,grid2y,iInterpEigVecGrid2(:,:,i))
    title(sprintf('Interpolated fine grid %d-st KL',i))
    
    subplot(1,3,3)
    surf(grid2x,grid2y,klGrid2(:,:,i))
    title(sprintf('Fine grid KL %d-st KL',i))
    
end

%% plot energy
figure(21)
plot(diag(eigVal1)*4)
hold on 
plot(diag(eigVal2))
hold off

%% Plot samples in line 
figure(3)
plot(Y2,'-')
hold on 
plot(interpY2,'--')
hold off

%% visualizing nKl effect
iSample=1;

for i=1:1:nKl

    iY2      =exp(klBasis2(:,1:i)*Theata(1:i,iSample)+interpMuXGrid2(:));
    iInterpY2=exp(interpKlBasis2(:,1:i)*Theata(1:i,iSample)+interpMuXGrid2(:));
    
    colorLimit=[min(iY2(:)) max(iY2(:))];
    
    figure(4)
    subplot(1,2,1)
    pcolor(grid2x,grid2y,reshape(iY2,size(grid2x)))
%     shading interp;
%     colormap jet;
    colorbar
    caxis(colorLimit);
    title(sprintf('Field sample nKL=%d',i))
    
    subplot(1,2,2)
    pcolor(grid2x,grid2y,reshape(iInterpY2,size(grid2x)))
%     shading interp;
%     colormap jet;
    colorbar
    caxis(colorLimit);
    title(sprintf('Interpolated Field sample nKL=%d',i))
    
    
    
%     contourf(X,Y,K(:,:,i))
% %     colormap(hot)
%     shading interp;
%     colorbar
%     title(sprintf('Permeability field'))
    frame(i)=getframe;
end

% [Phi,Lamda]=Kl()




