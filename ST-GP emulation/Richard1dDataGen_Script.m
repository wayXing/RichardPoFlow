% function [] = Richard1dDataGen_Script()
% 1d Spatial-temporial data generation using 1D Richard1d solver funtion.
% with heterogeneous permeability field as inputs.
%
% Input parameters:
%
% Output parameters:
%
% Examples:
%
% See also: 
% Author:   Wei Xing
% History:  02/08/2017 create
%
clear
close all

%% number of inputs samples and approximation accuracy
nSample=50;
nKl=10;


%% Solver Setup
% Spatial setup
lengthZ=100;
deltaZ=1;
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
lengthcale=lengthZ/10;
muY=0.0094; 
DeviationRatio=0.2;     %set DeviationRatio=10 to see dramatic results.
%     nKL=100;
klEnergyKeep=0.90;


[Z] = ndgrid(0:deltaZ:lengthZ);
%calculate distance matrix
distance = pdist(Z);
distanceMatrix = squareform(distance);
SigmaY=exp(-distanceMatrix./lengthcale) .*(muY*DeviationRatio)^2;      

% Conver to X covariance matrix and mean
SigmaX=log(SigmaY./(muY*muY')+ 1);
muX=log(muY)-diag(SigmaX)./2;

% KL/POD on X
[klBasis,klEigenValue,~] = svds(SigmaX,nKl);  % KL decomposition on covariance matrix via SVD/eigen decomposition

KlEnergy=diag(klEigenValue);
cumulatedKlEnergy= cumsum(KlEnergy)./sum(KlEnergy);
[~,miniNKl]=min(abs(cumulatedKlEnergy-klEnergyKeep))

% nKl=3;
sample= randn(nKl,nSample);

% Ks =exp(klBasis*sqrt(klEigenValue)*sample+repmat(muX,1,nSample));
Ksr=exp(klBasis(:,1:nKl)*sqrt(klEigenValue(1:nKl,1:nKl))*sample(1:nKl,:)+repmat(muX,1,nSample));
    

%% FOM on Kr
h=waitbar(0,'FOM on Ksr on progress');
for i=1:nSample
    mesh.Ks=Ksr(:,i);
    
    tic
    [hRecord(:,:,i),iteration(:,i)] = Richard1dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
    fomTimeCost(i)=toc;  
    nIterationCost(i)=sum(iteration(:,i));
    
    waitbar(i/nSample)
end



%% Plot
% figure(1)
% plot(cumulatedKlEnergy)
% hline =line([0,nKl],[klEnergyKeep,klEnergyKeep]);
% hline.Color = 'r';
% vline =line([nKl,nKl],[0,klEnergyKeep]);
% vline.Color = 'r';
% title(sprintf('Accumulated energy ration and truncation'))
% 
% 
% % subplot(2,2,2)
% figure(2)
% plot(nIterationCost)
% title(sprintf('number of iteration for each sample(on x axis)'))
% 
% 
% figure(3)
% plot(Ks)
% hold on 
% plot(Ksr)
% hold off
% title(sprintf('permeability field'))
% legend('All KL basis','Truncation KL basis')



%% Save
save('data','sample','Ksr','hRecord');
% load('data')




%% Auxiliary function
% function [n]=energy2n(allEigenvalues,energy)
% % KlEnergy=diag(klEigenValue);
% cumulatedKlEnergy= cumsum(allEigenvalues)./sum(energy);
% [~,n]=min(abs(cumulatedKlEnergy-klEnergyKeep));
% end

