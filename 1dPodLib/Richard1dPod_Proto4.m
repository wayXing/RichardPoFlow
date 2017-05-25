function [] = Richard1dPod_Proto4()
% Richars equation 1D pod solver testing file.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% Proto1: created from Richard1d_Demo3() show experiment 1d POD method.
% Proto2: take away function picardUpdate and use 1dPodLib
% Proto3: created from Proto3 to make DEIM POD
% Proto4: created from Proto3. This script runs true DEIM POD
%
% Input parameters:
%
% Output parameters:
%
% Examples:
%
% % Short description of example, followed by Matlab code line
% >> matlab code of example
%
% See also: 
% Author:   Wei Xing
% History:  22/05/2017  file created
%
clear
close all
%% Setup
% Spatial setup
lengthZ=100;
deltaZ=0.1;
nZ=lengthZ/deltaZ+1;

% Temporal setup
lengthTime=100;
deltaT=1;
nTime=lengthTime/deltaT;

% Iteration solver setup
nMaxIteration=50;
maxIteError=1;

%
podEnergy=0.995;

%nDeim
% nDeimK=30;
% nDeimC=30;

% nDeimK=nPOD;    %number of Deim basis for k term
% nDeimC=nPOD;    %number of Deim basis for c term


%% Initialize mesh
% [X,Z] = meshgrid(0:deltaX:lengthX,0:deltaZ:lengthZ);
[Z] = ndgrid(0:deltaZ:lengthZ);
mesh.lengthZ=lengthZ;
mesh.deltaZ=deltaZ;
mesh.nZ=nZ;


%%  Permeability field
scale=0.05;        % overall magnitude of the permeability field. decide the changing speed.
scale=0.0094; %recommand  
lengthcale=1;     %larger number means less stochastic (more correlation as one zooms in the 
                    %field) field. Thus gives smoother result.
                    
mesh.Ks=permeabilityField([Z(:)],lengthcale)*scale;
% mesh.Ks=ones(length(Z),1)*scale;

%% initial conditions and boundary value (DBC)
h_init=ones(nZ,1)*-61.5; %value for all initial points
h_init(1,1)=-20.7;       %value for top DBC
h_init(end,1)=-61.5;     %value for bottom DBC

mesh.dbcFlag=zeros(nZ,1);     %specify DBC location
mesh.dbcFlag(1)=1;
mesh.dbcFlag(end)=1;


%% Auxiliary variable   
% P=diag(dbcFlag);  picking up matrix
dbcIndex=find(mesh.dbcFlag);     %specify DBC index for later fitting in value
nodeIndex=find(~mesh.dbcFlag);   %specify free node index

nNode=sum(~mesh.dbcFlag);        %number of free node


%% FOM
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
%         mesh.K=kFunc(mesh.H);
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
    hRecord(:,t)=mesh.H;
    
end
fomTimeCostFom=toc  
nIterationFom

%% POD
% podEnergy=0.999;
hSnapShot=hRecord(:,:);

disp('POD decomposition process...')
[U,S,V]=svd(hSnapShot,'econ');

% [U,S,V]=svd(H);
energy=diag(S);
cumulatedEnergy= cumsum(energy)./sum(energy);
[~,nPOD]=min(abs(cumulatedEnergy-podEnergy))

% nPOD=80;
% U=U*sqrt(S);
V=U(:,1:nPOD);  %call the pod basis V


%% Naive ROM
mesh.H=h_init;
mesh.C=theataDifFunc(mesh.H);
% mesh.K=kFunc(mesh.H);
mesh.K=kFieldFunc(mesh.H,mesh.Ks);

nodeIndex=find(~mesh.dbcFlag);   %specify free node index   

nIterationPod=0;
tic
for t=1:nTime
    
    previousH=mesh.H;
    % Picard iteration
    for k=1:nMaxIteration 
        H0=mesh.H;  %preserved for iteration compare
        
        %update mesh value 
        mesh.C=theataDifFunc(mesh.H);
%         mesh.K=kFunc(mesh.H);
        mesh.K=kFieldFunc(mesh.H,mesh.Ks);
        
        [A,B]=picardAxbForm(mesh,previousH,deltaT);
        
        %POD decompose
        Ar=V(nodeIndex,:)'*A*V(nodeIndex,:);
        Br=V(nodeIndex,:)'*B;
        
        %solve linear equation
        hr=Ar\(Br);
        
        %reassemble
        h=V(nodeIndex,:)*hr;

        
        %update mesh value       
        mesh.H(nodeIndex)=h;
         
        %stopping criteria
        sseIte=sum((mesh.H(:)-H0(:)).^2);
        if sqrt(sseIte)<maxIteError 
            break 
        else 
            nIterationPod=nIterationPod+1;
        end
    end 
    hRecord2(:,t)=mesh.H;
    
end
podTimeCostFom=toc  
nIterationPod

%% DEIM nonlinear function 
% nDeimK=100;
% nDeimC=100;

nDeimK=nPOD;    %number of Deim basis for k term
nDeimC=nPOD;    %number of Deim basis for c term

%k
disp('DEIM decomposition for k...')
for t=1:size(hSnapShot,2)
    kRecord(:,t)=kFieldFunc(hSnapShot(:,t),mesh.Ks);
end
% [Vk,~,~]=svd(kRecord,'econ');
[Vk,~,~]=svds(kRecord,nDeimK);

[~,~,Pk] = DEIM(Vk);
Pk=Pk(:,1:nDeimK);
Vk=Vk(:,1:nDeimK);
Dk=Vk*inv(Pk'*Vk);  %DEIM basis

%c
disp('DEIM decomposition for c...')
cRecord=theataDifFunc(hSnapShot);
% [Vc,~,~]=svd(cRecord,'econ');
[Vc,~,~]=svds(cRecord,nDeimC);

[~,~,Pc] = DEIM(Vc);
Pc=Pc(:,1:nDeimC);
Vc=Vc(:,1:nDeimC);
Dc=Vc*inv(Pc'*Vc);  %DEIM basis


%% DEIM POD
%initilize system.
mesh.H=h_init;
mesh.C=theataDifFunc(mesh.H);
% mesh.K=kFunc(mesh.H);
mesh.K=kFieldFunc(mesh.H,mesh.Ks);

nodeIndex=find(~mesh.dbcFlag);   %specify free node index   

%Initialize ROM
[romMesh]=picardAxbRomInit(mesh,V,Dk,Dc);

VTV=V'*V;
Zh=V'*mesh.H;

nIterationDeimPod=0;
tic
for t=1:nTime   
    previousZh=Zh;
    % Picard iteration 
    for k=1:nMaxIteration 
        Zh0=Zh;  %preserved for iteration compare
         
        %update mesh value        
        Zk=kFieldFunc(Pk'*V*Zh,Pk'*mesh.Ks);
        Zc=theataDifFunc(Pc'*V*Zh);
           
        [Ar,Br]=picardAxbRomForm(romMesh,deltaT,previousZh,Zh,Zk,Zc);
             
        %solve linear equation
        Zh=Ar\(Br);
              
        %stopping criteria
%         mesh.H=V*Zh;
%         sseIte=sum((mesh.H(:)-H0(:)).^2);
        sseIte=(Zh-Zh0)'*VTV*(Zh-Zh0); 
        if sqrt(sseIte)<maxIteError 
            break 
        else 
            nIterationDeimPod=nIterationDeimPod+1;
        end
%         disp('.')
    end 
    ZhRecord(:,t)=Zh;
    
end
deimPodTimeCostFom=toc  

%after process
hRecord3=V*ZhRecord;
nIterationDeimPod



%% Plot
figure(1)
for t=1:1:nTime
    plot(hRecord(:,t))
    hold on 
    plot(hRecord2(:,t))
    plot(hRecord3(:,t))
    hold off
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;
end

% figure(2)
% for t=1:1:nTime
%     plot(V'*hRecord(:,t))
%     hold on 
%     plot(ZhRecord(:,t))
%     hold off
%     title(sprintf('ROM coifficients compare time=%i',t))
%     drawnow
%     frame(t)=getframe;
% end


  
end
%%



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

function result=kFunc(h)
rho=1.175e6;
r=4.74;
k_s=0.00944;

result=k_s.*rho./(rho+abs(h).^r);
end

function result=kFieldFunc(h,ks)
% h and k must be the same sizes
rho=1.175e6;
r=4.74;

result=ks.*rho./(rho+abs(h).^r);
end

%% permeability field generator
function Ks=permeabilityField(X,lengthcale)
% permeability generator (log-normal) given X coordinates and measure lengthscale.
%
% log(ks)=z~N(0,cov(x1,x2))
% cov[z_{12}]=cov(x1,x2)=exp(-|x1-x2|/c) and E[z]=0
%
% pointCoordinate=[X(:),Z(:)];

%larger number means less stochastic field. Thus less smooth.
% lengthcale=10; 

[nX,dimX]=size(X);

%calculate distance matrix
distance = pdist(X);
distanceMatrix = squareform(distance);

%calculate covariance matrix
covMatrix=exp(-distanceMatrix./lengthcale);    

% KL decomposition on covariance matrix via SVD/eigen decomposition
% [klBasis,klEigenValue] = eigs(covMatrix,nY*nX); 
[klBasis,klEigenValue,~] = svds(covMatrix,nX); 


% [nKlBasis,~]=sizes(klBasis);


%Generate independent normal samples 
seed=103;
rng(seed);
sample= randn(nX,1);

%make multivariate Gaussian distributions with samples. zero mean.
%Covariance specified though KL basis.
Ks=klBasis*sqrt(klEigenValue)*sample;
% Ks=reshape(Ks,nY,nX);

% a log (multi) normal permeability field
Ks=exp(Ks);

end

