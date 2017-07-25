function [] = Richard1dUQ_Proto2()
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
lengthZ=100;
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
lengthcale=lengthZ/10; %larger number means less stochastic (more correlation as one zooms in the 
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
klEnergyKeep=0.95;

sample= randn(nX,1);        %Random cofficient Sampling. Also the input in this case.
% [U,S,V]=svd(H);
KlEnergy=diag(klEigenValue);
cumulatedKlEnergy= cumsum(KlEnergy)./sum(KlEnergy);
[~,nKl]=min(abs(cumulatedKlEnergy-klEnergyKeep))

Ks =exp(klBasis*sqrt(klEigenValue)*sample).*scale;
Ksr=exp(klBasis(:,1:nKl)*sqrt(klEigenValue(1:nKl,1:nKl))*sample(1:nKl,1)).*scale;


%% FOM on K
mesh.Ks=Ks;

tic
[hRecord1,iteration1] = Richard1dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
fomTimeCostFom1=toc  
nIterationFom1=sum(iteration1)

%% FOM on Kr
mesh.Ks=Ksr;

tic
[hRecord2,iteration2] = Richard1dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
fomTimeCostFom2=toc  
nIterationFom2=sum(iteration2)


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
plot(iteration1)
hold on 
plot(iteration2)
hold off
title(sprintf('number of iteration at each time step'))
% legend('Full Kl','Truncated Kl')
legend(sprintf('Full Kl total=%i',sum(iteration1)),sprintf('Truncated Kl total=%i',sum(iteration2)))




figure(3)
plot(Ks)
hold on 
plot(Ksr)
hold off
title(sprintf('permeability field'))
legend('All KL basis','Truncation KL basis')

figure(4)
for t=1:1:nTime
    figure(4)
    plot(hRecord1(:,t))
    hold on 
    plot(hRecord2(:,t))
    hold off
    title(sprintf('time=%i',t))
%     legend('All KL basis','Truncation KL basis')
    drawnow
%     frame(t)=getframe;    %comment to save cpu time and memory
end




end





%% Auxiliary function
function [n]=energy2n(allEigenvalues,energy)
% KlEnergy=diag(klEigenValue);
cumulatedKlEnergy= cumsum(allEigenvalues)./sum(energy);
[~,n]=min(abs(cumulatedKlEnergy-klEnergyKeep));
end







