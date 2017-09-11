% function [] = stPcaEmulation()
% 1d Spatial-temporial data generation using 1D Richard1d solver funtion.
% with heterogeneous permeability field as inputs.
% load data by Richard1dDataGen_Script() and run PCA+ temporal emulation
%
% Input parameters:
%
% Output parameters:
%
% Examples:
%
% See also: 
% Author:   Wei Xing
% History:  02/09/2017 create
clear
close all
%%
filename='/Users/weix/Desktop/stgp-Data/data';
load(filename)

%% Define Train & Test data



%% PCA
EnergyKeepRatio=0.95;
[basis,eigValue,~]=svd(reshape(hRecord,[nZ,nTime*nSample]),'econ');

Energy=diag(eigValue);
cumulatedKlEnergyRatio= cumsum(Energy)./sum(Energy);
[~,nPca]=min(abs(cumulatedKlEnergyRatio-EnergyKeepRatio));

basis=basis(:,1:nPca);


%% Get latent z
% Z=pagefun(@mtimes,basis,hRecord);   %CUDA accerate

for i =1:nSample
    zRecord(:,:,i)=basis'*hRecord(:,:,i); 
    
    hRecordStar(:,:,i)=basis*zRecord(:,:,i);
end


%% Indepent GP 
for i=1:nPca
    
    
    
    
end





%% Compare Plot
ifPlot=1;
if ifPlot==1

nZShow=100;
zShow=1:round(nZ/nZShow):nZ;
for t=1:1:nTime
    figure(1)
    plot(squeeze( hRecord(zShow,t,:)),'-')
    ylim([-80,20])
    hold on 
    plot(squeeze( hRecordStar(zShow,t,:)),'--')
    hold off
    
    title(sprintf('time=%i',t))
%     legend('All KL basis','Truncation KL basis')
    drawnow
%     frame(t)=getframe;    %comment to save cpu time and memory
end

end

















%% Auxiliary function
% function [n]=energy2n(allEigenvalues,energy)
% KlEnergy=diag(klEigenValue);
% cumulatedKlEnergy= cumsum(allEigenvalues)./sum(energy);
% [~,n]=min(abs(cumulatedKlEnergy-klEnergyKeep));
% end