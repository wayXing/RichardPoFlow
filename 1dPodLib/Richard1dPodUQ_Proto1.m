function [] = Richard1dPodUQ_Proto1()
% UQ for Richars equation with random input 
%
% Richars equation 1D pod solver testing file.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% Proto1: introduce deim pod solver as a function. Also using
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

disp('KL decomposition on process...')
[nX,dimX]=size(Z);
nX=50;
[klBasis,klEigenValue,~] = svds(covMatrix,nX);  % KL decomposition on covariance matrix via SVD/eigen decomposition
% Vkl=klBasis*sqrt(klEigenValue);

%a log (multi) normal permeability field
%     Ks=exp(klBasis*sqrt(klEigenValue)*sample);

% randomCoief= randn(nX,1);

%% Make permeability field
%     klEnergyKeep=0.95;
% 
%     sample= randn(nX,1);        %Random cofficient Sampling. Also the input in this case.
%     % [U,S,V]=svd(H);
%     KlEnergy=diag(klEigenValue);
%     cumulatedKlEnergy= cumsum(KlEnergy)./sum(KlEnergy);
%     [~,nKl]=min(abs(cumulatedKlEnergy-klEnergyKeep))
% 
%     Ks =exp(klBasis*sqrt(klEigenValue)*sample).*scale;
%     Ksr=exp(klBasis(:,1:nKl)*sqrt(klEigenValue(1:nKl,1:nKl))*sample(1:nKl,1)).*scale;

sample= randn(size(klBasis,2),1);        %Random cofficient Sampling. Also the input in this case.
Ks =exp(klBasis*sqrt(klEigenValue)*sample).*scale;


%% FOM on K
mesh.Ks=Ks;

tic
[H_fom,iteration1] = Richard1dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
fomTimeCostFom1=toc  
nIterationFom1=sum(iteration1)


%% Prepare stage 
% POD
podEnergyKeep=0.995;
hSnapShot=H_fom(:,:);   %decide snapshot

disp('POD decomposition process...')
[U,S,~]=svd(hSnapShot,'econ');

% [U,S,V]=svd(H);
energy=diag(S);
cumulatedPodEnergy= cumsum(energy)./sum(energy);
[~,nPod]=min(abs(cumulatedPodEnergy-podEnergyKeep))

% nPOD=80;
% U=U*sqrt(S);
V=U(:,1:nPod);  %call the pod basis V

% DEIM nonlinear function 

% nDeimK=100;
% nDeimC=100;
nDeimK=nPod;    %number of Deim basis for k term
nDeimC=nPod;    %number of Deim basis for c term

%k
disp('DEIM decomposition for k...')
for t=1:size(hSnapShot,2)
    kRecord(:,t)=K(hSnapShot(:,t),Ks);
end
% [Vk,~,~]=svd(kRecord,'econ');
[Vk,~,~]=svds(kRecord,nDeimK);

[~,~,Pk] = DEIM(Vk);
Pk=Pk(:,1:nDeimK);
Vk=Vk(:,1:nDeimK);
Dk=Vk*inv(Pk'*Vk);  %DEIM basis

%c
disp('DEIM decomposition for c...')
cRecord=theataDif(hSnapShot);
% [Vc,~,~]=svd(cRecord,'econ');
[Vc,~,~]=svds(cRecord,nDeimC);

[~,~,Pc] = DEIM(Vc);
Pc=Pc(:,1:nDeimC);
Vc=Vc(:,1:nDeimC);
Dc=Vc*inv(Pc'*Vc);  %DEIM basis

%% Initilize ROM
[romMesh]=picardAxbRomInit(mesh,V,Dk,Dc);

tic
[H_pod,iteration2] = Richard1dPicardPodSolver(romMesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K,Pk,Pc);
podTimeCostFom1=toc  
nIterationFom1=sum(iteration1)


%% Plot
figure(1)

plot(cumulatedPodEnergy)
hline =line([0,nPod],[podEnergyKeep,podEnergyKeep]);
hline.Color = 'r';
vline =line([nPod,nPod],[0,podEnergyKeep]);
vline.Color = 'r';
title(sprintf('Accumulated energy ration and truncation for POD'))


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
for t=1:1:nTime
    figure(3)
    plot(H_fom(:,t))
    hold on 
    plot(H_pod(:,t))
    hold off
    title(sprintf('time=%i',t))
%     legend('All KL basis','Truncation KL basis')
    drawnow
%     frame(t)=getframe;    %comment to save cpu time and memory
end




end
