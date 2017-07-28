function [] = RIchard1dPod_Demo()
% Richard1d pod reduced order model solver demo funtion, 
% Uisng Pircards iteration on Dirichlet Boundary condition. The permeability field is
% heterogeneous.
%
% Input parameters:
%
% Output parameters:
%
% Examples:
%
% See also: 
% Author:   Wei Xing
% History:  27/07/2017  Document and modification
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
lengthcale=lengthZ/10;
muY=0.0094; 
DeviationRatio=0.4;     %set DeviationRatio=10 to see dramatic results.
nKl=50;
% klEnergyKeep=0.90;


[Z] = ndgrid(0:deltaZ:lengthZ);
%calculate distance matrix
distance = pdist(Z);
distanceMatrix = squareform(distance);
SigmaY=exp(-distanceMatrix./lengthcale) .*(muY*DeviationRatio)^2;      

% Conver to X covariance matrix and mean
SigmaX=log(SigmaY./(muY*muY')+ 1);
muX=log(muY)-diag(SigmaX)./2;

% KL/POD on X
[klBasis,klEigenValue,~] = svds(SigmaX,nZ);  % KL decomposition on covariance matrix via SVD/eigen decomposition

% KlEnergy=diag(klEigenValue);
% cumulatedKlEnergy= cumsum(KlEnergy)./sum(KlEnergy);
% [~,nKl]=min(abs(cumulatedKlEnergy-klEnergyKeep))

% nKl=3;
sample= randn(nZ,1);
Ks=exp(klBasis(:,1:nKl)*sqrt(klEigenValue(1:nKl,1:nKl))*sample(1:nKl,1)+muX);
    

%% FOM on K
mesh.Ks=Ks;

tic
[H_fom,iteration1] = Richard1dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
fomTimeCostFom1=toc  
nIterationFom1=sum(iteration1)


%% Prepare stage 
% POD
podEnergyKeep=0.995;
% podEnergyKeep=0.9999;

hSnapShot=H_fom(:,:);   %decide snapshot

disp('POD decomposition process...')
[U,S,~]=svd(hSnapShot,'econ');

% [U,S,V]=svd(H);
energy=diag(S);
cumulatedPodEnergy= cumsum(energy)./sum(energy);
[~,nPod]=min(abs(cumulatedPodEnergy-podEnergyKeep))

% nPod=40;
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
VdK=Vk*inv(Pk'*Vk);  %DEIM basis

%c
disp('DEIM decomposition for c...')
cRecord=theataDif(hSnapShot);
% [Vc,~,~]=svd(cRecord,'econ');
[Vc,~,~]=svds(cRecord,nDeimC);

[~,~,Pc] = DEIM(Vc);
Pc=Pc(:,1:nDeimC);
Vc=Vc(:,1:nDeimC);
VdC=Vc*inv(Pc'*Vc);  %DEIM basis

%% Deim POD
% Initilize ROM
[romMesh]=picardAxbRomInit(mesh,V,VdK,Pk,VdC,Pc);

tic
[H_pod,iteration2] = Richard1dPicardPodSolver(romMesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
podTimeCostFom1=toc  
nIterationFom1=sum(iteration1)


%% Plot
figure(1)
plot(Ks)
title(sprintf('Permeability field'))


figure(2)
plot(cumulatedPodEnergy)
hline =line([0,nPod],[podEnergyKeep,podEnergyKeep]);
hline.Color = 'r';
vline =line([nPod,nPod],[0,podEnergyKeep]);
vline.Color = 'r';
title(sprintf('Accumulated energy ration and truncation for POD'))


% subplot(2,2,2)
figure(3)
plot(iteration1)
hold on 
plot(iteration2)
hold off
title(sprintf('number of iteration at each time step'))
% legend('Full Kl','Truncated Kl')
legend(sprintf('Full Kl total=%i',sum(iteration1)),sprintf('Truncated Kl total=%i',sum(iteration2)))


figure(4)
for t=1:1:nTime
    figure(4)
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
