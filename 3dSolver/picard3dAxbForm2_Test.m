function [] = picard3dAxbForm2_Test()
% picard3dAxbForm2 Test function
% The function focus on fix Dirichlet BC.
%
% Input parameters:
%
% Output parameters:
%
% Examples:
%
% See also: 
% Author:   Wei Xing
% History:  31/05/2017  file created
%           25/07/2017  Document and modification
% 
% Log:
% Version1.0 -the process need further vectorizition and parellel to speed
%             up and save memory.
%            -The way boundary condition are stored and indicated need
%            improvements.
%            -the way coordinate are defined as z,x,y NOT conventional x,y,z (because z
%            is more of interest for Richard's equation) Thus conventional 
%            plotting becomes difficult.

%% Setup
% Spatial setup
% Spatial setup
lengthZ=40;
deltaZ=4;
nZ=lengthZ/deltaZ+1;

lengthX=40;
deltaX=4;
nX=lengthX/deltaX+1;

lengthY=40;
deltaY=4;
nY=lengthY/deltaY+1;

% Temporal setup
lengthT=100;
deltaT=1;
nTime=lengthT/deltaT;

% Iteration solver setup
nMaxIteration=50;
maxIteError=1;

[Z,X,Y] = ndgrid(0:deltaZ:lengthZ,0:deltaX:lengthX,0:deltaY:lengthY);


%% Define Boundary conditions
% Option 1       
% Make nodeIndex matrix. It indicate the sequence nodes are accessed and the type of node.  
%     nodeIndex=zeros(nZ,nX);
%     nodeIndex(2:end-1,2:end-1)=reshape(uint32(1:(nZ-2)*(nX-2)), (nZ-2), (nX-2));    
% 
%     nodeInFieldIndex=find(nodeIndex);
% 
%     %%% 
%     %initial state
%     H_init=zeros(nZ,nX);
%     H_init(2:end-1,2:end-1)=-61.5;
% 
%     %BC
%     % bcLeft=ones(nNodeZ,1)*-20.7;
%     % bcRight=ones(nNodeZ,1)*-61.5;
%     % bcTop=ones(nNodeX,1)*-20.7;
%     % bcBottom=ones(nNodeX,1)*-61.5;
% 
% 
%         % A more interesting setup
%         bcLeft=ones(nZ,1)*-20.7;
%         bcRight=ones(nZ,1)*-20.7;
%         bcTop=ones(nX,1)*-20.7;
%         bcBottom=ones(nX,1)*-24.7;
% 
% 
%     H_init(1,1:end)=bcTop;
%     H_init(end,1:end)=bcBottom;
% 
%     H_init(1:end,1)=bcLeft;
%     H_init(1:end,end)=bcRight;

% Option 2 Neumann BC on sides and Dirichlet conditions on Top and bottom
    nodeIndex=zeros(nZ,nX,nY);
    nodeIndex(2:end-1,1:end,1:end)=reshape(uint32(1:(nZ-2)*nX*nY), (nZ-2),nX,nY);  
    
    nodeInFieldIndex=find(nodeIndex);

    nodeIndex(:,:,1)=-nodeIndex(:,:,1);
    nodeIndex(:,:,end)=-nodeIndex(:,:,end);
    
    nodeIndex(:,1,2:end-1)=-nodeIndex(:,1,2:end-1);
    nodeIndex(:,end,2:end-1)=-nodeIndex(:,end,2:end-1);
    
    
%% Define initial condition
    H_init=ones(nZ,nX,nY)*-61.5;
    H_init(1,:,:)=ones(nX,nY)*-20.7;
    H_init(end,:,:)=ones(nX,nY)*-61.5;
    
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


%% Define Permeability field input 
scale=0.05;        % overall magnitude of the permeability field. decide the changing speed.
lengthScale=100;      %larger number means less stochastic (more correlation as one zooms in the 
                    %field) field. Thus gives smoother result.
nKl=50;             %number of KL baisi/cofficients

[nAllNode,nDimension]=size([Z(:),X(:),Y(:)]);
% nKl=nAllNode;

%calculate distance matrix
distance = pdist([Z(:),X(:),Y(:)]);
distanceMatrix = squareform(distance);

%calculate covariance matrix
covMatrix=exp(-distanceMatrix./lengthScale);    

% KL decomposition on covariance matrix via SVD/eigen decomposition
[klBasis,klEigenValue,~] = svds(covMatrix,nKl); 


% Generate independent normal samples 
seed=100;
rng(seed);
sample= randn(nKl,1);

%make multivariate Gaussian distributions with samples. zero mean.
%Covariance specified though KL basis.
% Ks=klBasis*sqrt(klEigenValue)*sample;
% Ks=reshape(Ks,nY,nX);
% a log (multi) normal permeability field
% Ks=exp(Ks);

Ks =exp(klBasis*sqrt(klEigenValue)*sample).*scale;

Ks=reshape(Ks,nZ,nX,nY);

%Plot 
bubbleScale=100;
scatter3(X(:),Y(:),Z(:),Ks(:)*bubbleScale,Ks(:)*bubbleScale)


%         %Homogenous Ks
%         muY=0.0094; 
%         Ks=repmat(muY,nZ,nX,nY); 

%% Initilize mesh 
mesh.deltaZ=deltaZ;
mesh.nZ=nZ;
mesh.deltaX=deltaX;
mesh.nX=nX;
mesh.deltaY=deltaY;
mesh.nY=nY;

mesh.nodeIndex=nodeIndex;
mesh.nNode=length(nodeIndex(nodeIndex~=0));

mesh.H=H_init;      %impose initial condition.

%% MAIN
% nNode=length(nodeIndex(nodeIndex~=0));
% H=H_init;       %impose initial condition.
h=waitbar(0,'On progress');
tic
for t=1:nTime
    H_PreviousTime= mesh.H;
    for k=1:nMaxIteration       
        H0=mesh.H;
        
        %update mesh value 
        mesh.C=theataDif(mesh.H);
        mesh.K=K(mesh.H,Ks);
        
%         [A,B] = PicardFdm(H_PreviousTime);
        [A,B]=picard3dAxbForm2(mesh,H_PreviousTime,deltaT);

        hFree=A\(-B);
                
%         H(find(nodeIndex))=hFree;       %pay extra attention to ordering
        mesh.H(nodeInFieldIndex)=hFree; 
%         mesh.H=H;
        
        sseIte=sum((mesh.H(:)-H0(:)).^2);
        if sqrt(sseIte)<maxIteError 
            break 
        end
        
    end
    
    hRecord(:,:,:,t)=mesh.H;
    waitbar(t/nTime)
end
toc
close(h)

%% Plotting 
% H=mesh.H;

% figure(1)
% scatter3(X(:),Y(:),Z(:),Ks(:)*10000,Ks(:)*10000)

% end time pressure
figure(2)
scatter3(X(:),Y(:),Z(:),abs(mesh.H(:)),abs(mesh.H(:)))
title('end time pressure')

% figure(5)
% p = patch(isosurface(H,-60));
% isonormals(H,p)
% p.FaceColor = 'red';
% p.EdgeColor = 'none';
% daspect([1 1 1])
% view(3); 
% axis([nX,nY,nZ ])

% preasure head front surface propogation
figure(3)
for t=1:nTime
    clf
    Ht=abs(hRecord(:,:,:,t));
    p=patch(isosurface(Ht,40));
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    
    p=patch(isosurface(Ht,50));
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';  

    p=patch(isosurface(Ht,60));
    p.FaceColor = 'blue';
    p.EdgeColor = 'none';  
    
    view(-60,40)
    title(sprintf('pressure front propogation. t=%i',t))
    axis([1,nX,1,nY,1,nZ])
    
    camlight 
    lighting gouraud

    
    drawnow
    frame(t)=getframe;
end

xslice = [nX]; 
yslice = [nY]; 
zslice = [0:nZ/2:nZ];
figure(3)
slice(mesh.H,nX,1:5:nZ,1)


% preasure head propogation in square volume
figure(4)
for t=1:nTime
%     surf(X(:,:,sliceY),Z(:,:,sliceY),TheataRecord(:,:,sliceY,t))
%     shading interp;
%     slice(TheataRecord(:,:,:,t),nX,nY,0:nZ/4:nZ)
    slice(hRecord(:,:,:,t),nX-3,1:5:nZ,1)
    view(-60,40)
    title(sprintf('pressure. t=%i',t))
    drawnow
    frame(t)=getframe;
end

% preasure head propogation in sphere volume
figure(5)
for t=1:4:nTime
    hVector=(hRecord(:,:,:,t)+100)*20;    
%     scatter3(X(:),Y(:),Z(:),(hVector(:)),(Ks(:)*10000),'fill')
    scatter3(X(:),Y(:),Z(:),(hVector(:)),(hVector(:)),'fill')
%     shading interp;
    title(sprintf('pressure. t=%i',t))
    drawnow
    frame(t)=getframe;
    
end




