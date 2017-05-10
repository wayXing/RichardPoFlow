function [] = Richard1d_Demo3()
% Richars equation 1D solver tester.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% Demo3: use shirft matrix operation to formulate function of ax=b;
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
% History:  06/05/2017  file created
%
%% Parameter Initial 
%     % Set standard initialization parameter
%        InitOptStandard = [0.2, 0.25, 0];
%        % InitRand = 0.3 (use normal random inoculation with 30% of domain of variables),
%        % InitNindKeep = 0.2 (keep 20% maximal)
%        % InitNindUniform = 0 (create no random individuals)
%     % Check input parameters
%        % At least the first 2 input parameters are necessary
%        if nargin < 2,  error('Not enough input parameters (at least 2 parameters - individuals and Nind)'); end
%        if isnan(PopInit), PopInit = []; end
%        % The 3. input parameter is optional, default is []
%        if nargin < 3,  VLUB = []; end
%     % the 4. input parameter is optional and can contain multiple options
%        if nargin < 4, InitOpt = []; end
%        if isnan(InitOpt), InitOpt = []; end
%        % When too many options are contained in InitOpt, issue a warning and 
%        % shorten the parameter vector
%        if length(InitOpt) > length(InitOptStandard), 
%           InitOpt = InitOpt(1:length(InitOptStandard)); 
%           warning(' Too many parameters in InitOpt! InitOpt was shortened.');
%        end
tic
%% Setup
% Spatial setup
lengthZ=10000;
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
% mesh.Ks=permeabilityField([Z(:)],lengthcale)*scale;


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


%% Main 
mesh.H=h_init;

mesh.C=theataDifFunc(mesh.H);
mesh.K=kFunc(mesh.H);


for t=1:nTime
    
    mesh=picardUpdate(mesh,deltaT,nMaxIteration,maxIteError);   
    TheataRecord(:,:,t)=mesh.H;
    
end

toc  
figure(1)
for t=1:nTime
    plot(TheataRecord(:,:,t))
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;

end

  
end
%%
function mesh=picardUpdate(mesh,deltaT,nMaxIteration,maxIteError)
%update mesh value using Picards iteration
  
    previousH=mesh.H;
    
    %main
    
    for k=1:nMaxIteration 
        H0=mesh.H;  %preserved for iteration compare
        
        %update mesh value 
        mesh.C=theataDifFunc(mesh.H);
        mesh.K=kFunc(mesh.H);
        
        [A,B]=picardAxbForm2(mesh,previousH,deltaT);
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


end

function [A,B]=picardAxbForm(mesh,previousH,deltaT)
%calculate linear system equation Ax=B of free node.
% This function is probabilly the fastest as not auxiliary shirft
% matrix/Picking up matrix is formed explicitly. All done by matrix
% indexing and circshirft

    % Auxiliary variable   
    dbcIndex=find(mesh.dbcFlag);     %specify DBC index for later fitting in value
    nodeIndex=find(~mesh.dbcFlag);   %specify free node index
    nNode=sum(~mesh.dbcFlag);        %number of free node   

    deltaZ=mesh.deltaZ;
    nZ=mesh.nZ;

    %main
%     C=theataDifFunc(mesh.H);
%     K=kFunc(mesh.H);
    
    C=mesh.C;
    K=mesh.K;
    
%     K=kFieldFunc(mesh.H,mesh.Ks);

    centerDiag = (2.*K+ circshift(K,1)+circshift(K,-1))/(2*deltaZ^2)+C/deltaT;  %A_center diagonal %first and last elements are meaningless                    
    upDiag     = (K+ circshift(K,1))/(-2*deltaZ^2);                                %A_up     diagonal %first and last elements are meaningless                  
    downDiag   = (K+ circshift(K,-1))/(-2*deltaZ^2);                               %A_down   diagonal %first and last elements are meaningless 

    B          = -(circshift(K,-1)-circshift(K,1))/(2*deltaZ)+previousH.*C/deltaT;

    %make spare A                    
    %Todo this part may be improved using better sparse diag
    Amethod=1;
    switch Amethod 
        case 1      %Very fast
            A_all=spdiags(centerDiag,0,nZ,nZ) +circshift(spdiags(upDiag,0,nZ,nZ),[0,-1]) +circshift(spdiags(downDiag,0,nZ,nZ),[0,1]);                                          
        case 2
            A_all=sparse (diag(centerDiag,0) +circshift(diag(upDiag,0),[0,-1]) +circshift(diag(downDiag,0),[0,1]));       
        otherwise 
    end

    %pick free node and componsate for dbc involved 

    B=B(nodeIndex)-A_all(nodeIndex,dbcIndex)*mesh.H(dbcIndex);
    A=A_all(nodeIndex,nodeIndex);

end


% POD picardAxbForm
function [A,B]=picardAxbForm2(mesh,previousH,deltaT)
%calculate linear system equation Ax=B of free node with shirft matrix

    % Auxiliary variable   
    dbcIndex=find(mesh.dbcFlag);     %specify DBC index for later fitting in value
    nodeIndex=find(~mesh.dbcFlag);   %specify free node index
    nNode=sum(~mesh.dbcFlag);        %number of free node   

    deltaZ=mesh.deltaZ;
    nZ=mesh.nZ;
   
    C=mesh.C;
    K=mesh.K;
    
    % P=diag(dbcFlag);
    % P=spdiags(dbcFlag,0,nZ,nZ)
    
    % Calculation
    
    Amethod=1;
    switch Amethod
        case 1  %vary fast
            UpShift1Eye =circshift(speye(nZ),[-1,0]);
            lowShift1Eye=circshift(speye(nZ),[1,0]);
        case 2  %fast              
            UpShift1Eye =circshift(spdiags(ones(nZ,1),0,nZ,nZ),[-1,0]);
            lowShift1Eye=circshift(spdiags(ones(nZ,1),0,nZ,nZ),[1,0]);
        case 3  %DONT use 
        % If sparse matrix is not adopted. later calculation would be slow!
%             UpShirft1Eye=spdiags(ones(nZ),0,nZ,nZ);
            UpShift1Eye=diag(ones(nZ,1));                  %spdiags and diag takes similar time but different memory 
            UpShift1Eye=circshift(UpShift1Eye,[-1,0]);    %circshift is quick and use almost no time 
            lowShift1Eye=diag(ones(nZ,1));
            lowShift1Eye=circshift(lowShift1Eye,[1,0]);
    end
    
%     K=kFieldFunc(mesh.H,mesh.Ks);
    
    Amethod=1;
    switch Amethod
        case 1  
            %write 3 diagonal band  %very fast. as avoid matrix 
            %calculation and use a element wise operation.
            centerDiag = (2.*K+ lowShift1Eye*K+UpShift1Eye*K)/(2*deltaZ^2)+C/deltaT;       %A_center diagonal %first and last elements are meaningless                    
            upDiag     = (K+ lowShift1Eye*K)/(-2*deltaZ^2);                                %A_up     diagonal %first and last elements are meaningless                  
            downDiag   = (K+ UpShift1Eye*K)/(-2*deltaZ^2);                                 %A_down   diagonal %first and last elements are meaningless         
        case 2     
            %rewrite 3 diagonal band  %fast. do matrix calculation though
            %using sparse operations.
            centerDiag = (2.*speye(nZ)+lowShift1Eye+UpShift1Eye)./(2*deltaZ^2)*K+C/deltaT;   
            upDiag     = (speye(nZ)+lowShift1Eye)./(2*deltaZ^2)*K;
            downDiag   = (speye(nZ)+UpShift1Eye)./(2*deltaZ^2)*K;
    end
    
    B          = -(UpShift1Eye*K-lowShift1Eye*K)/(2*deltaZ)+previousH.*C/deltaT;

    %make spare metrix A                    
    %Todo this part may be improved using better sparse diag
    Amethod=4;
    switch Amethod 
        case 1      %Very fast
            A_all=spdiags(centerDiag,0,nZ,nZ) +circshift(spdiags(upDiag,0,nZ,nZ),[0,-1]) +circshift(spdiags(downDiag,0,nZ,nZ),[0,1]);                                          
        case 2      %Slow. Always NOT creat full matrix and then sparse it.
            A_all=sparse (diag(centerDiag,0) +circshift(diag(upDiag,0),[0,-1]) +circshift(diag(downDiag,0),[0,1]));   
        case 3      %not necessary but show a more clear formulation
            left1Eye=circshift(spdiags(ones(nZ),0,nZ,nZ),[0,-1]);
            right1Eye=circshift(spdiags(ones(nZ),0,nZ,nZ),[0,1]);
            A_all=spdiags(centerDiag,0,nZ,nZ) +spdiags(upDiag,0,nZ,nZ)*left1Eye + spdiags(downDiag,0,nZ,nZ)*right1Eye;   
        case 4
            A_all=spdiags(centerDiag,0,nZ,nZ) +spdiags(upDiag,0,nZ,nZ)*lowShift1Eye + spdiags(downDiag,0,nZ,nZ)*UpShift1Eye;  
        case 5
%             A_all=spdiags(centerDiag,0,nZ,nZ) +spdiags(upDiag,0,nZ,nZ)*lowShift1Eye + spdiags(downDiag,0,nZ,nZ)*UpShift1Eye;  
        case 6
            A_all=centerDiag*speye(nZ) +upDiag*speye(nZ)*lowShift1Eye + downDiag*speye(nZ)*lowShift1Eye;  
        otherwise 
    end

    %pick free node and componsate for dbc involved 

    B=B(nodeIndex)-A_all(nodeIndex,dbcIndex)*mesh.H(dbcIndex);
    A=A_all(nodeIndex,nodeIndex);

end




%% 
function mesh=picardUpdate2(mesh,deltaT,nMaxIteration,maxIteError)
% update mesh value using Picards iteration including calculate linear system equation Ax=B of free node.

    % Auxiliary variable   
    dbcIndex=find(mesh.dbcFlag);     %specify DBC index for later fitting in value
    nodeIndex=find(~mesh.dbcFlag);   %specify free node index
    nNode=sum(~mesh.dbcFlag);        %number of free node   
    
    deltaZ=mesh.deltaZ;
    nZ=mesh.nZ;
    
    previousH=mesh.H;
    
    %main
    for k=1:nMaxIteration 
    
        H0=mesh.H;  %preserved for iteration compare
                
        C=theataDifFunc(mesh.H);
        K=kFunc(mesh.H);
%         K=kFieldFunc(mesh.H,mesh.Ks);

        centerDiag = (2.*K+ circshift(K,1)+circshift(K,-1))/(2*deltaZ^2)+C/deltaT;  %A_center diagonal %first and last elements are meaningless                    
        upDiag     = (K+ circshift(K,1))/(-2*deltaZ^2);                                %A_up     diagonal %first and last elements are meaningless                  
        downDiag   = (K+ circshift(K,-1))/(-2*deltaZ^2);                               %A_down   diagonal %first and last elements are meaningless 

        B      = (circshift(K,-1)-circshift(K,1))/(2*deltaZ)-previousH.*C/deltaT;

        %make spare A                    
        %Todo this part maybe improved using better sparse diag
        Amethod=1;
        switch Amethod 
            case 1      %Very fast
                A_all=spdiags(centerDiag,0,nZ,nZ) +circshift(spdiags(upDiag,0,nZ,nZ),[0,-1]) +circshift(spdiags(downDiag,0,nZ,nZ),[0,1]);                                          
            case 2
                A_all=sparse (diag(centerDiag,0) +circshift(diag(upDiag,0),[0,-1]) +circshift(diag(downDiag,0),[0,1]));       
            otherwise 
        end

        %pick free node and componsate for dbc involved 

        B=B(nodeIndex)+A_all(nodeIndex,dbcIndex)*mesh.H(dbcIndex);
        A=A_all(nodeIndex,nodeIndex);


        %solve linear equation
        h=A\(-B);
        mesh.H(nodeIndex)=h;
        
        %stopping criteria
        sseIte=sum((mesh.H(:)-H0(:)).^2);
        if sqrt(sseIte)<maxIteError 
            break 
        end
        
    end
       
            
end



%% 
function picardAbForm()
% different way to update mesh value using Picards iteration including calculate linear system equation Ax=B of free node.

            imethod=3;
            switch imethod 
                case 1      %simple, easy to understand but slow program
                    
                    indexMatrix=zeros(nZ,1);
                    indexMatrix(2:end-1)=1:nZ-2;                
                    i=1;
                    
                    for j =2:nZ-1

                        kHalfUp   =(K(j,i)+K(j-1,i))/2;
                        kHalfDown =(K(j,i)+K(j+1,i))/2;

                        cCenter=C(j,i);

                        wUp   = -kHalfUp  ./deltaZ^2;
                        wDown = -kHalfDown./deltaZ^2;

                        wCenter=cCenter/deltaTime-wUp-wDown;
                        b=(kHalfDown-kHalfUp)/deltaZ-H_PreviousTime(j,i)*cCenter/deltaTime;

                        indexCenter=indexMatrix(j,i);
                        indexUp=indexMatrix(j-1,i);
                        indexDown=indexMatrix(j+1,i);

                        %Check BC and modify 
                        %check if top BC
                        if j==2 
                            indexUp=0;
                            b=b+wUp*H(j-1,i);
                        end

                        %check if bottom BC
                        if j==nZ-1 
                            indexDown=0;
                            b=b+wDown*H(j+1,i);
                        end

                        if indexUp>0 A(indexCenter,indexUp)=wUp; end
                        if indexDown>0 A(indexCenter,indexDown)=wDown; end

                        A(indexCenter,indexCenter)=wCenter;
                        B(indexCenter,1)=b;

                    end
                    
                case 2         %% Method II semi-Vectorize code (prepare for Method III)        
                    
                    Atemp=zeros(nZ);
                    Btemp=zeros(nZ,1);
%                     nZ=nNode;

                    for i=2:nZ-1    % for all non-boundary points
                       Atemp(i,i)   = (2*K(i)+ K(i-1)+K(i+1))/(2*deltaZ^2)+C(i)/deltaTime;
                       Atemp(i,i-1) = (K(i)+ K(i-1))/(-2*deltaZ^2);
                       Atemp(i,i+1) = (K(i)+ K(i+1))/(-2*deltaZ^2); 
                       Btemp(i)     = (K(i+1)-K(i-1))/(2*deltaZ)-H_PreviousTime(i)*C(i)/deltaTime;

                    end

                    A=Atemp(nodeIndex,nodeIndex);
                    B=Btemp(nodeIndex)+Atemp(nodeIndex,dbcIndex)*H(dbcIndex);
                    
                case 3      %% Method III Vectorize code    
                    
                    centerDiag = (2.*K+ circshift(K,1)+circshift(K,-1))/(2*deltaZ^2)+C/deltaTime;  %A_center diagonal %first and last elements are meaningless                    
                    upDiag     = (K+ circshift(K,1))/(-2*deltaZ^2);                                %A_up     diagonal %first and last elements are meaningless                  
                    downDiag   = (K+ circshift(K,-1))/(-2*deltaZ^2);                               %A_down   diagonal %first and last elements are meaningless 
                        
                    Btemp      = (circshift(K,-1)-circshift(K,1))/(2*deltaZ)-H_PreviousTime.*C/deltaTime;
                                      
                    %make spare A                    
                    %Todo this part maybe improved using better sparse diag
                    Amethod=1;
                    switch Amethod 
                        case 1      %Very fast
                            Atemp=spdiags(centerDiag,0,nZ,nZ) +circshift(spdiags(upDiag,0,nZ,nZ),[0,-1]) +circshift(spdiags(downDiag,0,nZ,nZ),[0,1]);                                          
                        case 2
                            Atemp=sparse (diag(centerDiag,0) +circshift(diag(upDiag,0),[0,-1]) +circshift(diag(downDiag,0),[0,1]));       
                        otherwise 
                    end
                                                                        
                    %pick free node and componsate for dbc involved 
                    A=Atemp(nodeIndex,nodeIndex);
                    B=Btemp(nodeIndex)+Atemp(nodeIndex,dbcIndex)*H(dbcIndex);
             
                otherwise 
                    error('no such method')
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

