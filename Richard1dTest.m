function [] = Richard1dTest()
% Richars equation 1D solver tester.
% The function focus on fix Dirichlet BC.
% This function serves as a simple tester for all Richard solver developed in this
% project.
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
% History:  03/05/2017  file created
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
lengthZ=40;
deltaZ=1;
nZ=lengthZ/deltaZ+1;

% Temporal setup
lengthTime=300;
deltaTime=1;
nTime=lengthTime/deltaTime;

% Iteration solver setup
nMaxIteration=1000;
miniIteError=0.1;


%% Initialize mesh
% [X,Z] = meshgrid(0:deltaX:lengthX,0:deltaZ:lengthZ);
[Z] = ndgrid(0:deltaZ:lengthZ);

%%  Permeability field
lengthcale=100;     
%larger number means less stochastic (more correlation as one zooms in the 
%field) field. Thus gives smoother result.
scale=0.005;
% overall magnitude of the permeability field. decide the changing speed.

Ks=permeabilityField([Z(:)],lengthcale)*scale;


%% initial conditions and boundary value (DBC)
h_init=ones(nZ,1)*-61.5;               %value for all initial points

h_init(1,1)=-20.7;       %value for top DBC
h_init(end,1)=-61.5;     %value for bottom DBC


dbcFlag=zeros(nZ,1);
dbcFlag(1)=1;
dbcFlag(end)=1;

%Auxiliary variable
% P=diag(dbcFlag);

dbcIndex=find(dbcFlag);
nodeIndex=find(~dbcFlag);

nNode=sum(~dbcFlag);


%% Main

imethod=3;

H=h_init;
for t=1:nTime
    
    H_PreviousTime=H;
    
    for k=1:nMaxIteration 
    
        H0=H;
        
        C=theataDifFunc(H);
        K=kFunc(H);
        
            %% Method I
            %initialize
            A=zeros(nNode);
            B=zeros(nNode,1);

            %Assemble Ax+B=0;

            switch imethod 
                case 1 
                    
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

%                         indexCenter=j;
%                         indexUp=j-1;
%                         indexDown=j+1;

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
                    
                case 2          %% Method II semi-Vectorize code                                                  
                    Atemp=zeros(nZ);
                    Btemp=zeros(nZ,1);
%                     nZ=nNode;

                    for i=2:nZ-1    % for all non-boundary points

                %        w_center(i)= (2*K(i)+ K(i-1)+K(i+1))/(4*deltaZ^2)+C(i)/deltaT;
                %        w_up(i)    = (K(i)+ K(i-1))/(-2*deltaZ^2);
                %        w_down(i)  = (K(i)+ K(i+1))/(-2*deltaZ^2); 
                %        b(i)       = (K(i+1)-K(i-1))/(2*deltaZ);

                       Atemp(i,i)   = (2*K(i)+ K(i-1)+K(i+1))/(2*deltaZ^2)+C(i)/deltaTime;
                       Atemp(i,i-1) = (K(i)+ K(i-1))/(-2*deltaZ^2);
                       Atemp(i,i+1) = (K(i)+ K(i+1))/(-2*deltaZ^2); 
                       Btemp(i)     = (K(i+1)-K(i-1))/(2*deltaZ)-H_PreviousTime(i)*C(i)/deltaTime;

                    end

%                     dbcFlag=zeros(nZ,1);
%                     dbcFlag(1)=1;
%                     dbcFlag(end)=1;
% 
%                     P=diag(dbcFlag);
% 
% 
%                     dbcIndex=find(dbcFlag);
%                     nodeIndex=find(~dbcFlag);

                    A=Atemp(nodeIndex,nodeIndex);
                    B=Btemp(nodeIndex)+Atemp(nodeIndex,dbcIndex)*H(dbcIndex);
                    
                    
                    
                    
                case 3      %% Method III Vectorize code    
                    
                    Ac = (2.*K+ circshift(K,1)+circshift(K,-1))/(2*deltaZ^2)+C/deltaTime;   %first and last elements are meaningless                    
%                         Ac([1,end])=[0;0];      %make explicit for first and last elements to be meaningless. This line is NOT necessary for the code to run.   

                    Au = (K+ circshift(K,1))/(-2*deltaZ^2);                         
%                         Au([end,end-1])=[0;0];      %make explicit for first and last elements to be meaningless. This line is NOT necessary for the code to run.       
                    
                    Ad = (K+ circshift(K,-1))/(-2*deltaZ^2);                            
%                         Ad([1,2])=[0;0];      %make explicit for first and last elements to be meaningless. This line is NOT necessary for the code to run. 
                        
                    Btemp  = (circshift(K,-1)-circshift(K,1))/(2*deltaZ)-H_PreviousTime.*C/deltaTime;
                                                                       
                    %make spare A                    
                    %Todo this part should be improve by directly creat
                    %sparse band matrix.
                    
                    Amethod=1;
                    switch Amethod 
                        case 1      %Very fast
                            Atemp=spdiags(Ac,0,nZ,nZ) +circshift(spdiags(Au,0,nZ,nZ),[0,-1]) +circshift(spdiags(Ad,0,nZ,nZ),[0,1]);                                          
                        case 2
                            Atemp=sparse (diag(Ac,0) +circshift(diag(Au,0),[0,-1]) +circshift(diag(Ad,0),[0,1]));       
                        otherwise 
%                       % Try(1)    Fail
%                         Atemp=spdiags(Ac,0) +circshift(spdiags(Au,0),[0,-1]) +circshift(spdiags(Ad,0),[0,1]);
%                       % Try(2)    Fail  
%                         Atemp=speye(nZ)*Ac +circshift(speye(nZ)*Au,[0,-1]) +circshift(speye(nZ)*Ad,[0,1]);
                    end
                                                
%                     diagAu=circshift(diag(Au,0),[0,-1]);
             
                        %Todo this part should be improve by directly creat
                        %sparse band matrix.
                        %real sparse

%                       % Try(3)    Work & time saving  
                        

                  
                    %select free node
                    A=Atemp(nodeIndex,nodeIndex);
                    B=Btemp(nodeIndex)+Atemp(nodeIndex,dbcIndex)*H(dbcIndex);
                    
                    
                otherwise 
                    error('no such method')
            end
 

            
            
%         %trime for BC
%         indexIfKnown=~any(A);
%         indexIfUnknown=any(A);
% 
%         A=A(indexIfUnknown,indexIfUnknown);
%         B=B(indexIfUnknown);


        h=A\(-B);
%         H(2:end-1)=h;
        H(nodeIndex)=h;
    
        
        sseIte=sum((H(:)-H0(:)).^2);
        if sqrt(sseIte)<miniIteError 
            break 
        end
        
    end
    TheataRecord(:,:,t)=H;
    
end
    
    toc
    
    figure(1)
    % plot(H_init);

    for t=1:nTime
        plot(TheataRecord(:,:,t))
        title(sprintf('time=%i',t))
        drawnow
        frame(t)=getframe;

    end

  
end


%% 
function D=vec2dig(v,D,i)
% add vector to a (sparse) diagal matrix at i diagonal 
% size(v,1) must equal to size(D,1). 

    Dc=diag(v);
    




end





function [A2,b2]=Picard_v1(H,H_PreviousTime,deltaT,K,C,nZ)

    A=zeros(nZ);
    b=zeros(nZ,1);
    
    for i=2:nZ-1    % for all non-boundary points
        
%        w_center(i)= (2*K(i)+ K(i-1)+K(i+1))/(4*deltaZ^2)+C(i)/deltaT;
%        w_up(i)    = (K(i)+ K(i-1))/(-2*deltaZ^2);
%        w_down(i)  = (K(i)+ K(i+1))/(-2*deltaZ^2); 
%        b(i)       = (K(i+1)-K(i-1))/(2*deltaZ);
       
       A(i,i)  = (2*K(i)+ K(i-1)+K(i+1))/(4*deltaZ^2)+C(i)/deltaT;
       A(i,i-1)  = (K(i)+ K(i-1))/(-2*deltaZ^2);
       A(i,i+1)  = (K(i)+ K(i+1))/(-2*deltaZ^2); 
       b(i)       = (K(i+1)-K(i-1))/(2*deltaZ);
       
    end

    dbcFlag=zeros(nZ);
    dbcFlag(1)=1;
    dbcFlag(end)=1;
    
    P=diag(dbcFlag);
    
    
    dbcIndex=find(dbcFlag);
    nodeIndex=~dbcIndex;
    
    A2=A(nodeIndex,nodeIndex);
    b2=b(nodeIndex)-A(nodeIndex,dbcIndex)*H(dbcIndex);


end








function [A,B]=axbPicard(H,H_PreviousTime,K,C)
%% Assemble Ax+B=0 with New Vectorization 
% H,K,C are value at all grid point (in order to calculate coifficients)
% size (A) ~= size(H);
%         H0=H;
    H_all(2:end-1,2:end-1)=H;

    C_all=theataDifFunc(H_all);
    C=C_all(2:end-1,2:end-1);

    K_all=kFunc(H_all);
    K=K_all(2:end-1,2:end-1);

    zdiffK_all=diff(K_all,1,1);
    xdiffK_all=diff(K_all,1,2);

    wUp   = -1./deltaZ^2 .* (K - zdiffK_all(1:end-1,2:end-1)./2);
    wDown = -1./deltaZ^2 .* (K + zdiffK_all(2:end,  2:end-1)./2);

    wLeft = -1./deltaX^2 .* (K - xdiffK_all(2:end-1,1:end-1)./2);
    wRight= -1./deltaX^2 .* (K + xdiffK_all(2:end-1,2:end)./2);

    wCenter=C./deltaTime-wUp-wDown-wLeft-wRight;

    b= (zdiffK_all(2:end, 2:end-1)+ zdiffK_all(1:end-1,2:end-1))./2 ./deltaZ ...
       -H_PreviousTime .* C ./ deltaTime;

    %update BC neighbour influence to Ax+B=0              
    B= b(:) + wUp(:)   .*topDbcValue...
            + wDown(:) .*bottomDbcValue...
            + wLeft(:) .*leftDbcValue...
            + wRight(:).*rightDbcValue;    

    wUp    = wUp(:)   .* (topDbcValue==0); 
    wDown  = wDown(:) .* (bottomDbcValue==0); 
    wLeft  = wLeft(:) .* (leftDbcValue==0); 
    wRight = wRight(:).* (rightDbcValue==0); 
    wCenter=wCenter(:);

    %Way I. heavy time cost. Not sparse
%         A=  diag(wUp(2:end),1)+diag(wDown(1:end-1),-1)...
%             +diag(wLeft(nNodeZ+1:end),-nNodeZ)+diag(wRight(1:end-nNodeZ),nNodeZ)...       
%             +diag(wCenter,0);

    %Way II.
    band=zeros(nNodeZ*nNodeX,5);
    band(1:end-nNodeZ,1)= wLeft(nNodeZ+1:end);
    band(1:end-1,2)= wDown(1:end-1);
    band(:,3)= wCenter;
    band(1+1:end,4)= wUp(1+1:end);
    band(1+nNodeZ:end,5)= wRight(1:end-nNodeZ);

    A = spdiags(band,[-nNodeZ,-1,0,1,nNodeZ],nNodeZ*nNodeX,nNodeZ*nNodeX);

end
    
  
function [A,b]=RichardsPicard(h0,hPreviousTime,deltaT,K,C)

        H_all(2:end-1,2:end-1)=H; 
        C_all=theataDifFunc(H_all);
    
    
    
    
    
    
end









%%
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

