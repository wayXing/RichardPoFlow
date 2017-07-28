function [] = Richard1dPodUQ_Proto2()
% UQ for Richars equation with random input 
%
% Richars equation 1D pod solver testing file.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% Proto1: introduce deim pod solver as a function. Also using
%         inline function to describe non-linear term.
% Proto2: created from Proto1 and Richard1dUQ_Proto3. Batch UQ inout field
%         for 1d Pod with independent ideal basis.
%
% Input parameters:
%
% Output parameters:
%
% See also: 
%
% Author:   Wei Xing
% History:  29/05/2017  file created
%
% tips:     -Ks with higher randomness would lead to highly variate result
%           and thus unstable results for pod.
%           -Long time simulation sould lead to unstable results. (evolved
%            basis would be required/Very large number of basis)
%           -Long time pod simulation most fails at the begining. start with
%           Fom then switch to FOM would help stabalize.  
%           -Fine grid e.g. deltaZ=0.01 would show significant cpu saving.
%           -Fine grid also improve stability


%
clear
close all
%% Solver Setup
% Spatial setup
lengthZ=100;
deltaZ=0.1;
nZ=lengthZ/deltaZ+1;

% Temporal setup
lengthT=200;
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
%         % scale=0.05;  
%         nKl=5;
%         scale=0.0094;            %recommand value from paper 
%         lengthcale=0.75*lengthZ; %larger number means less stochastic (more correlation as one zooms in the 
%                                  %field) field. Thus gives smoother result.
% 
%         [Z] = ndgrid(0:deltaZ:lengthZ);
%         %calculate distance matrix
%         distance = pdist(Z);
%         distanceMatrix = squareform(distance);
% 
%         covMatrix=exp(-distanceMatrix./lengthcale);    %calculate covariance matrix
% 
%         % [nX,dimX]=size(Z);
%         % nKl=nX;
% 
%         disp('KL decomposition for Ks...')
%         [klBasis,klEigenValue,~] = svds(covMatrix,nKl);  % KL decomposition on covariance matrix via SVD/eigen decomposition
%         % Vkl=klBasis*sqrt(klEigenValue);
% 
%         %a log (multi) normal permeability field
%         %     Ks=exp(klBasis*sqrt(klEigenValue)*sample);
% 
%         % randomCoief= randn(nX,1);
% 
%          %% Make permeability field
%         nSample=10;
%         % klEnergyKeep=0.95;
% 
%         sample= randn(nKl,nSample);        %Random cofficient Sampling. Also the input in this case.
%         % [U,S,V]=svd(H);
%         % KlEnergy=diag(klEigenValue);
%         % cumulatedKlEnergy= cumsum(KlEnergy)./sum(KlEnergy);
%         % [~,nKl]=min(abs(cumulatedKlEnergy-klEnergyKeep))
% 
%         Ks =exp(klBasis*sqrt(klEigenValue)*sample).*scale;
%         % Ksr=exp(klBasis(:,1:nKl)*sqrt(klEigenValue(1:nKl,1:nKl))*sample(1:nKl,:)).*scale;
        
%% Define and Decompose the permeability input field II
lengthScale=lengthZ*0.2; %larger number means less stochastic (more correlation as one zooms in the 
nKl=30;
nSample=10;

% KsMean=0.0094;
% KsVar= (KsMean*0.3)^2;

GaussianMean= log(0.0094);
GaussianVar = (GaussianMean*0.2)^2;       


[Z] = ndgrid(0:deltaZ:lengthZ);
%calculate distance matrix
distance = pdist(Z);
distanceMatrix = squareform(distance);

covMatrix=exp(-distanceMatrix./lengthScale);    %calculate correlation matrix 
covMatrix=GaussianVar*covMatrix;                %calculate covariance  matrix 

[klBasis,klEigenValue,~] = svds(covMatrix,nKl);  % KL decomposition on covariance matrix via SVD/eigen decomposition

% Make permeability field
% sample= randn( size(klBasis,2),1);               %Sampling from a normal distribution
sample= randn(nKl,nSample);                              %Sampling from a normal distribution

Ks = (klBasis*sqrt(klEigenValue)*sample)+GaussianMean;  %Multivariate Gaussian
Ks = exp(Ks);                                           %Multivariate  Log-normal

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

%% Prepare stage 
% POD
% podEnergyKeep=0.995;
nPod=40;
   
h=waitbar(0,'independent Pod and Deimon K&C on progress');
for i=1:nSample    
    hSnapShot=H_uq1(:,:,i);   %decide snapshot

   
%     [U,S,~]=svd(hSnapShot,'econ');
%     energy=diag(S);
%     cumulatedPodEnergy= cumsum(energy)./sum(energy);
%     [~,nPod]=min(abs(cumulatedPodEnergy-podEnergyKeep))
%     % U=U*sqrt(S);
%     V(:,:,i)=U(:,1:nPod);  %call the pod basis V
    
    [V_uq(:,:,i),S,~]=svds(hSnapShot,nPod);
    

    % DEIM nonlinear function 

    % nDeimK=100;
    % nDeimC=100;
    nDeimK=nPod;    %number of Deim basis for k term
    nDeimC=nPod;    %number of Deim basis for c term

    %k
%     disp('DEIM decomposition for k...')
    for t=1:size(hSnapShot,2)
        kRecord(:,t)=K(hSnapShot(:,t),Ks(:,i));
    end
    % [Vk,~,~]=svd(kRecord,'econ');
    [Vk,~,~]=svds(kRecord,nDeimK);

    [~,~,Pk] = DEIM(Vk);
    Pk=Pk(:,1:nDeimK);
    Vk=Vk(:,1:nDeimK);
    Dk_uq(:,:,i)=Vk*inv(Pk'*Vk);  %DEIM basis
    
    %special treatment to store sparse logical matrix Pk
    [iRow,iColume]=find(Pk);
    PkiRow_uq(:,i)=iRow;
%     Pk=sparse(iRow,iColume,ones(nDeimK,1),nZ,nDeimK);     %recovery

    %c
%     disp('DEIM decomposition for c...')
    cRecord=theataDif(hSnapShot);
    % [Vc,~,~]=svd(cRecord,'econ');
    [Vc,~,~]=svds(cRecord,nDeimC);

    [~,~,Pc] = DEIM(Vc);
    Pc=Pc(:,1:nDeimC);
    Vc=Vc(:,1:nDeimC);
    Dc_uq(:,:,i)=Vc*inv(Pc'*Vc);  %DEIM basis

    %special treatment to store sparse logical matrix Pk
    [iRow,iColume]=find(Pc);
    PciRow_uq(:,i)=iRow;
%     Pc=sparse(iRow,iColume,ones(nDeimC,1),nZ,nDeimC);     %recovery
    
    waitbar(i/nSample)
end
close(h)

%% Deim POD
h=waitbar(0,'Deim pod on Ks on progress');


for i=1:nSample
    
    % Initilize ROM
    Pk=sparse(PkiRow_uq(:,i),iColume,ones(nDeimK,1),nZ,nDeimK);     %recovery form storage
    Pc=sparse(PciRow_uq(:,i),iColume,ones(nDeimC,1),nZ,nDeimC);     %recovery form storage
    mesh.Ks=Ks(:,i);
    
%     mesh.H=H_uq1(:,1,i); % use fom to start
    tic 
    [romMesh]=picardAxbRomInit(mesh,V_uq(:,:,i),Dk_uq(:,:,i),Pk,Dc_uq(:,:,i),Pc);
    
    tCostPodInit(i,1)=toc;
    
    tic
    [H_pod,iteration2] = Richard1dPicardPodSolver(romMesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
    
%     tCost2Total(i,1)=toc;
    tCost2(i,1)=toc;
    
    iTera2(i,2)=sum(iteration2);
    H_uq2(:,:,i)=H_pod;
    
    waitbar(i/nSample)
end
close(h)

% tCost2=tCost2Total-tCostPodInit;
tCost2Total=tCost2+tCostPodInit;

sum(tCost1)
sum(tCost2)

%% UQ process
mu_H_uq1 =mean(H_uq1,3);
var_H_uq1=std(H_uq1,0,3);
mid_H_uq1=median(H_uq1,3);

mu_H_uq2 =mean(H_uq2,3);
var_H_uq2=std(H_uq2,0,3);
mid_H_uq2=median(H_uq2,3);


%% plot

%% show basis
v1=reshape(V_uq(:,1,:),nZ,nSample);
figure(1)
plot(v1)

%show cpu time
figure(1)
plot([tCost1,tCost2,tCost2Total]);
title(sprintf('cpu time cost'))
legend('fom','pod run time','pod total time (with model building time')



nZShow=100;
zShow=1:round(nZ/nZShow):nZ;
figure(3)
for t=1:1:nTime
    figure(3)
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



figure(2)

nZShow=100;
zShow=1:round(nZ/nZShow):nZ;
figure(2)
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



end
