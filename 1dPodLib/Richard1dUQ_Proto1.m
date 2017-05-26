function [] = Richard1dUQ_Proto1()
% UQ for Richars equation with random input 
%
% Richars equation 1D pod solver testing file.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% Proto1: Runs FOM only. It compares results to show if the KL-decomposition
% works well. only one sample.
%
% Input parameters:
%
% Output parameters:
%
% See also: 
%
% Author:   Wei Xing
% History:  24/05/2017  file created
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


% Iteration solver setup
nMaxIteration=50;
maxIteError=1;

% 
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


%% Define and Decompose the permeability field
% scale=0.05;  
scale=0.0094; %recommand  
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
sample= randn(nX,1);        %Random cofficient Sampling. Also the input in this case.

%% Make permeability field
klEnergyKeep=0.95;

% [U,S,V]=svd(H);
KlEnergy=diag(klEigenValue);
cumulatedKlEnergy= cumsum(KlEnergy)./sum(KlEnergy);
[~,nKl]=min(abs(cumulatedKlEnergy-klEnergyKeep))

Ks =exp(klBasis*sqrt(klEigenValue)*sample).*scale;
Ksr=exp(klBasis(:,1:nKl)*sqrt(klEigenValue(1:nKl,1:nKl))*sample(1:nKl,1)).*scale;



%% Define non-linear function 
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

rho=1.175e6;
r=4.74;
% kFromhKs @(h) Ks.*rho./(rho+abs(h).^r);
K = @(h) Ks.*rho./(rho+abs(h).^r);

theata    = @(h)  alpha.*(theata_s-theata_r)/(alpha+abs(h).^beta)+theata_r;
theataDif = @(h) -alpha.*(theata_s-theata_r).*-1.*(alpha+abs(h).^beta).^(-2).*abs(h).^(beta-1);

%% FOM on K
mesh.Ks=Ks;

mesh.H=h_init;
mesh.C=theataDifFunc(mesh.H);
mesh.K=kFieldFunc(mesh.H,mesh.Ks);

hSnapShot=[];
nIterationFom=0;
tic

tStep=1:deltaT:lengthT;

for iT=1:length(tStep)
    
    previousH=mesh.H;
      % Picard iteration
    for k=1:nMaxIteration 
        H0=mesh.H;  %preserved for iteration compare
        
        %update mesh value 
%         mesh.C=theataDifFunc(mesh.H);
% %         mesh.K=kFunc(mesh.H);
%         mesh.K=kFieldFunc(mesh.H,mesh.Ks);
        
        mesh.C=theataDif(mesh.H);
        mesh.K=K(mesh.H);
        
        
        [A,B]=picardAxbForm(mesh,previousH,deltaT);
        %solve linear equation
        h=A\(B);
        
        %update mesh value
        nodeIndex=find(~mesh.dbcFlag);   %specify free node index   
        mesh.H(nodeIndex)=h;
         
        hSnapShot=[hSnapShot,mesh.H];
        %stopping criteria
        sseIte=sum((mesh.H(:)-H0(:)).^2);
        if sqrt(sseIte)<maxIteError 
            iterationRecordFom(iT)=k;
            break 
        else 
            nIterationFom=nIterationFom+1;
        end
    end 
    hRecord1(:,iT)=mesh.H;
    
end
fomTimeCostFom=toc  
nIterationFom



%% FOM on Kr
mesh.Ks=Ksr;

mesh.H=h_init;
mesh.C=theataDifFunc(mesh.H);
mesh.K=kFieldFunc(mesh.H,mesh.Ks);

hSnapShot=[];
nIterationFom=0;
tic
for t=1:nTime
    
    previousH=mesh.H;
      % Picard iteration
    for k=1:nMaxIteration 
        H0=mesh.H;  %preserved for iteration compare
        
        %update mesh value 
        mesh.C=theataDifFunc(mesh.H);
        mesh.K=kFieldFunc(mesh.H,mesh.Ks);
        
        [A,B]=picardAxbForm(mesh,previousH,deltaT);
        %solve linear equation
        h=A\(B);
        
        %update mesh value
        nodeIndex=find(~mesh.dbcFlag);   %specify free node index   
        mesh.H(nodeIndex)=h;
         
        hSnapShot=[hSnapShot,mesh.H];
        %stopping criteria
        sseIte=sum((mesh.H(:)-H0(:)).^2);
        if sqrt(sseIte)<maxIteError 
            iterationRecordFom(t)=k;
            break 
        else 
            nIterationFom=nIterationFom+1;
        end
    end 
    hRecord2(:,t)=mesh.H;
    
end
fomTimeCostFom=toc  
nIterationFom


%% Plot
figure(1)
plot(Ks)
hold on 
plot(Ksr)
hold off
title(sprintf('permeability field'))
legend('All KL basis','Truncation KL basis')



figure(2)
for t=1:1:nTime
    plot(hRecord1(:,t))
    hold on 
    plot(hRecord2(:,t))
    hold off
    title(sprintf('time=%i',t))
%     legend('All KL basis','Truncation KL basis')
    drawnow
%     frame(t)=getframe;
end



end


%% Other functions
function theata=theataFunc(h)
theataS=0.287;
theataR=0.075;
alpha=1.611e6;
beta=3.96;

result=alpha.*(theataS-theataR)/(alpha+abs(h).^beta)+theataR;
end

function theataDif=theataDifFunc(h)
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

theataDif=-alpha.*(theata_s-theata_r).*-1.*(alpha+abs(h).^beta).^(-2).*abs(h).^(beta-1);

end

% function result=kFunc(h)
% rho=1.175e6;
% r=4.74;
% k_s=0.00944;
% 
% result=k_s.*rho./(rho+abs(h).^r);
% end

function result=kFieldFunc(h,ks)
% h and k must be the same sizes
rho=1.175e6;
r=4.74;

result=ks.*rho./(rho+abs(h).^r);
end

%% Auxiliary function
function [n]=energy2n(allEigenvalues,energy)
% KlEnergy=diag(klEigenValue);
cumulatedKlEnergy= cumsum(allEigenvalues)./sum(energy);
[~,n]=min(abs(cumulatedKlEnergy-klEnergyKeep));
end







