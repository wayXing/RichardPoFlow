function [] = Richard1dUQ_Proto3()
% UQ for Richars equation with random input 
%
% Richars equation 1D pod solver testing file.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% Proto1: Runs FOM only. It compares results to show if the KL-decomposition
%         works well. only one sample.
% Proto2: created from Proto1. introduce solver as a function. Also using
%         inline function to describe non-linear term.
% Proto3: created from Proto2. Batch input to do UQ.
%
% Input parameters:
%
% Output parameters:
%
% See also: 
%
% Author:   Wei Xing
% History:  26/05/2017  file created
%
clear
close all
%% Solver Setup
% Spatial setup
lengthZ=10;
deltaZ=0.1;
nZ=lengthZ/deltaZ+1;

% Temporal setup
lengthT=300;
deltaT=1;
nTime=lengthT/deltaT;
% tStep=1:deltaT:lengthTime;


%Solver iteration setup
nMaxIteration=50;
maxIteError=1;

% update mesh structure
mesh.lengthZ=lengthZ;
mesh.deltaZ=deltaZ;
mesh.nZ=nZ;

%% initial conditions and boundary value (DBC)
h_init=ones(nZ,1)*-61.5; %value for all initial points
h_init(1,1)=-20.7;       %value for top DBC
h_init(end,1)=-61.5;     %value for bottom DBC
mesh.H=h_init;

mesh.dbcFlag=zeros(nZ,1);     %specify DBC location
mesh.dbcFlag(1)=1;
mesh.dbcFlag(end)=1;

%% Define for C and K non-linear function
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

rho=1.175e6;
r=4.74;
% kFromhKs @(h) Ks.*rho./(rho+abs(h).^r);
% K = @(h) Ks.*rho./(rho+abs(h).^r);
K = @(h,Ks) Ks.*rho./(rho+abs(h).^r);

theata    = @(h)  alpha.*(theata_s-theata_r)/(alpha+abs(h).^beta)+theata_r;
theataDif = @(h) -alpha.*(theata_s-theata_r).*-1.*(alpha+abs(h).^beta).^(-2).*abs(h).^(beta-1);

%% Define and Decompose the permeability input field
% scale=0.05;  
scale=0.0094;          %recommand value from paper 
lengthcale=lengthZ/4; %larger number means less stochastic (more correlation as one zooms in the 
                       %field) field. Thus gives smoother result.
              
[Z] = ndgrid(0:deltaZ:lengthZ);
%calculate distance matrix
distance = pdist(Z);
distanceMatrix = squareform(distance);

covMatrix=exp(-distanceMatrix./lengthcale);    %calculate covariance matrix

[nX,dimX]=size(Z);
[klBasis,klEigenValue,~] = svds(covMatrix,nX);  % KL decomposition on covariance matrix via SVD/eigen decomposition
% Vkl=klBasis*sqrt(klEigenValue);

%a log (multi) normal permeability field
%     Ks=exp(klBasis*sqrt(klEigenValue)*sample);

% randomCoief= randn(nX,1);

%% Make permeability field
nSample=10;
klEnergyKeep=0.95;

sample= randn(nX,nSample);        %Random cofficient Sampling. Also the input in this case.
% [U,S,V]=svd(H);
KlEnergy=diag(klEigenValue);
cumulatedKlEnergy= cumsum(KlEnergy)./sum(KlEnergy);
[~,nKl]=min(abs(cumulatedKlEnergy-klEnergyKeep))

Ks =exp(klBasis*sqrt(klEigenValue)*sample).*scale;
Ksr=exp(klBasis(:,1:nKl)*sqrt(klEigenValue(1:nKl,1:nKl))*sample(1:nKl,:)).*scale;

% scale=0.000094;
% Ks =exp(klBasis*sqrt(klEigenValue)*sample+scale);
% Ksr=exp(klBasis(:,1:nKl)*sqrt(klEigenValue(1:nKl,1:nKl))*sample(1:nKl,:)+scale);

%% FOM on K
% define non-linear term
h=waitbar(0,'FOM on Ks on progress');
for i=1:nSample
%     K = @(h) Ks(:,i).*rho./(rho+abs(h).^r);  %!CALL this function every time uisng new permeability field Ks
    mesh.Ks=Ks(:,i);
    tic
    [H,iteration1] = Richard1dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
    tCost1(i,1)=toc;
    
    iTera1(i,1)=sum(iteration1);
    H_uq1(:,:,i)=H;
    
    waitbar(i/nSample)
end
close(h)

%% FOM on Kr
% define non-linear term


h=waitbar(0,'FOM on Ksr on progress');
for i=1:nSample
%     K = @(h) Ksr(:,i).*rho./(rho+abs(h).^r);  %!CALL this function every time uisng new permeability field Ks/Ksr
    mesh.Ks=Ksr(:,i);
    tic
    [H,iteration2] = Richard1dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
    tCost2(i,1)=toc;
    
    iTera2(i,1)=sum(iteration2);
    H_uq2(:,:,i)=H;
    
    waitbar(i/nSample)
end
close(h)

%% UQ process
mu_H_uq1 =mean(H_uq1,3);
var_H_uq1=std(H_uq1,0,3);
mid_H_uq1=median(H_uq1,3);

mu_H_uq2 =mean(H_uq2,3);
var_H_uq2=std(H_uq2,0,3);
mid_H_uq2=median(H_uq2,3);

%% Plot
figure(1)

plot(cumulatedKlEnergy)
hline =line([0,nKl],[klEnergyKeep,klEnergyKeep]);
hline.Color = 'r';
vline =line([nKl,nKl],[0,klEnergyKeep]);
vline.Color = 'r';
title(sprintf('Accumulated energy ration and truncation'))


% subplot(2,2,2)
figure(2)
plot(iTera1)
hold on 
plot(iTera2)
hold off
title(sprintf('number of iteration at each sample run'))
% legend('Full Kl','Truncated Kl')
legend(sprintf('Full Kl total=%i',sum(iTera1)),sprintf('Truncated Kl total=%i',sum(iTera2)))


figure(3)
plot(Ks,'-')
hold on 
plot(Ksr,'--')
hold off
title(sprintf('permeability field'))
legend('All KL basis','Truncation KL basis')

zShow=1:1:nZ;
figure(4)
for t=1:1:nTime
    figure(4)
    errorbar(zShow,mu_H_uq1(zShow,t),var_H_uq1(zShow,t),'-')
    hold on
    plot(zShow,mid_H_uq1(zShow,t),'-')
    
    errorbar(zShow,mu_H_uq2(zShow,t),var_H_uq2(zShow,t),'--')
    plot(zShow,mid_H_uq2(zShow,t),'--')
    hold off
    title(sprintf('Mean Variance and Median @t=%i',t))
%     legend('All KL basis','Truncation KL basis')
    drawnow
%     frame(t)=getframe;    %comment to save cpu time and memory
end


zShow=1:1:nZ;
figure(5)
for t=1:1:nTime
    figure(5)
    plot(squeeze( H_uq1(zShow,t,:)),'-')
    hold on 
    plot(squeeze( H_uq2(zShow,t,:)),'--')
    hold off
    ylim([-80,20])
    
    title(sprintf('time=%i',t))
%     legend('All KL basis','Truncation KL basis')
    drawnow
%     frame(t)=getframe;    %comment to save cpu time and memory
end





% bar3(squeeze(H_uq1(:,3,:)))






sum(tCost1)
sum(tCost2)



end





%% Auxiliary function
function [n]=energy2n(allEigenvalues,energy)
% KlEnergy=diag(klEigenValue);
cumulatedKlEnergy= cumsum(allEigenvalues)./sum(energy);
[~,n]=min(abs(cumulatedKlEnergy-klEnergyKeep));
end







