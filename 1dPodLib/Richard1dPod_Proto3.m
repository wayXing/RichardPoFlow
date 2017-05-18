function [] = Richard1dPod_Proto3()
% Richars equation 1D pod solver testing file.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% Proto1: created from Richard1d_Demo3() show experiment 1d POD method.
% Proto2: take away function picardUpdate and use 1dPodLib
% Proto3: created from Proto3 to make DEIM POD
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
% History:  18/05/2017  file created
%
tic
%% Setup
% Spatial setup
lengthZ=100;
deltaZ=1;
nZ=lengthZ/deltaZ+1;

% Temporal setup
lengthTime=300;
deltaT=1;
nTime=lengthTime/deltaT;

% Iteration solver setup
nMaxIteration=1000;
maxIteError=0.1;


%% Initialize mesh
% [X,Z] = meshgrid(0:deltaX:lengthX,0:deltaZ:lengthZ);
[Z] = ndgrid(0:deltaZ:lengthZ);
mesh.lengthZ=lengthZ;
mesh.deltaZ=deltaZ;
mesh.nZ=nZ;


%%  Permeability field
scale=0.005;        % overall magnitude of the permeability field. decide the changing speed.
lengthcale=10;     %larger number means less stochastic (more correlation as one zooms in the 
                    %field) field. Thus gives smoother result.
mesh.Ks=permeabilityField([Z(:)],lengthcale)*scale;


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
         
        %stopping criteria
        sseIte=sum((mesh.H(:)-H0(:)).^2);
        if sqrt(sseIte)<maxIteError 
            break 
        end
    end 
    TheataRecord(:,t)=mesh.H;
    
end
fomTimeCostFom=toc  


%% POD
podEnergy=0.999;

[U,S,V]=svd(TheataRecord,'econ');
% [U,S,V]=svd(H);
energy=diag(S);
cumulatedEnergy= cumsum(energy)./sum(energy);
[~,nPOD]=min(abs(cumulatedEnergy-podEnergy))

V=U(:,1:nPOD);  %call the pod basis V


%% Naive ROM
mesh.H=h_init;
mesh.C=theataDifFunc(mesh.H);
% mesh.K=kFunc(mesh.H);
mesh.K=kFieldFunc(mesh.H,mesh.Ks);

nodeIndex=find(~mesh.dbcFlag);   %specify free node index   

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
        end
    end 
    TheataRecord2(:,t)=mesh.H;
    
end
podTimeCostFom=toc  

%% DEIM POD











%% Plot
figure(1)
for t=1:nTime
    plot(TheataRecord(:,t))
    hold on 
    plot(TheataRecord2(:,t))
    hold off
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;

end


  
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
seed=101;
rng(seed);
sample= randn(nX,1);

%make multivariate Gaussian distributions with samples. zero mean.
%Covariance specified though KL basis.
Ks=klBasis*sqrt(klEigenValue)*sample;
% Ks=reshape(Ks,nY,nX);

% a log (multi) normal permeability field
Ks=exp(Ks);

end

