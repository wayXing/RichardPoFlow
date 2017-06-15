function [] = Richard1dMlPodUqInit_Proto1()
% Initilize UQ for Richars equation with random input using mchaine learning & POD
% generate training testing dataset
%
% Richars equation 1D pod solver testing file.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% 
% Proto1: Created form Richard1dMlPodUQ_Proto1() in idPod
%
% Input parameters:
%
% Output parameters:
%
% See also: 
%
% Author:   Wei Xing
% History:  15/06/2017  file created
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
        
%% Define and Decompose the permeability input field II
lengthScale=lengthZ*0.2; %larger number means less stochastic (more correlation as one zooms in the 
nKl=30;
nSample=200;

% KsMean=0.0094;
% KsVar= (KsMean*0.3)^2;

GaussianMean= log(0.0014);
GaussianVar = (GaussianMean*0.4)^2;       


[Z] = ndgrid(0:deltaZ:lengthZ);
%calculate distance matrix
distance = pdist(Z);
distanceMatrix = squareform(distance);

covMatrix=exp(-distanceMatrix./lengthScale);    %calculate correlation matrix 
covMatrix=GaussianVar*covMatrix;                %calculate covariance  matrix 

[klBasis,klEigenValue,~] = svds(covMatrix,nKl);  % KL decomposition on covariance matrix via SVD/eigen decomposition

% Make permeability field
% sample= randn( size(klBasis,2),1);               %Sampling from a normal distribution
% sample= randn(nKl,nSample);                              %Sampling from a normal distribution
sample= lhsdesign(nSample,nKl)';                         %Sampling using latin hypercube

% sample(1:1,:)=0.5*ones(nSample,1);

Ks = (klBasis*sqrt(klEigenValue)*sample)+GaussianMean;  %Multivariate Gaussian
Ks = exp(Ks);                                          


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

%% Storage
filename='/Users/weix/Desktop/Data/data1';
% [status,msg] = mkdir('/uesrs/desktop/Data/data1.mat')
save(filename);







end
