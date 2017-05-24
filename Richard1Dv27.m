function [  ] = Richard1Dv27()
% 1D Richards equation with Dirichlet/Neumann boundary condition 
% h based Richards equation
%
% First edition: Weix 21/04/2017 
%
% point type 0: normal 
%  		  	 1: Dirichlet BC
%		  	 2: Neuman BC 
% 
% TODO BD inside. for now Neumann BC can only placed in the boundary
% update vesrion v2.3: to get rid of fdmGrid class.
%                      change the way points are calculated
% update vesrion v2.5: Updata Neumann boundary condition and the way nodes
%                      are accessed.(step 1. scan all inner point and indentify if they are 
%                      next to BC point. step 2.scan all outer boundary point)
% update vesrion v2.5.2: add support to permeability field.
% update vesrion v2.7: improve readibility and flexibility (sacrefice some performance)
%                      Compared to v2.5. all nodes are scaned at one time
%                      then identidy its type and neighbour.

tic
%% Grid structure Initial
lengthZ=1000;
deltaZ=1;

nZ=lengthZ/deltaZ+1;    %number of Z index total 

Z=[0:deltaZ:lengthZ]';

h=zeros(nZ,1);
nodeIndex=zeros(nZ,1);

Ks=permeabilityField(Z)*0.01;
%% costumize initial value and BC

%option 1 dbc
%     nodeIndex(2:end-1)=1:nZ-2;
%     h=ones(nZ,1).*-60;
% 
%     ifdbcNode(1)=1;
%     h(1)=-20;
% 
%     % ifdbcNode(end)=1;
%     h(end)=-60;
%     
%     ifNbc=logical(zeros(nZ,1));     %no Neumann BC

%option 2 dbc + nbc
    nodeIndex(2:end)=1:nZ-1;
    h=ones(nZ,1).*-60;
    
    ifdbcNode=logical(zeros(nZ,1)); 
    ifdbcNode(1)=1;
    h(1)=-20;

    ifNbc=logical(zeros(nZ,1));         %default zero flux
    ifNbc(end)=1;
    % nbcValue=2;       certain flux


%% solver setup
lengthTime=300;
deltaTime=1;
nTime=lengthTime/deltaTime+1;

% Iteration setup
nMaxIteration=1000;
miniIteError=0.01;


%% MAIN
nNode=length(nodeIndex(nodeIndex~=0));

for t=1:nTime
      
    h_PreviousTime= h;
    for k=1:nMaxIteration 
        
        h0=h;
        
        [A,B] = PicardFdm(h_PreviousTime);
        hFree=A\(-B);
                
        h(find(nodeIndex))=hFree;       %pay extra attention to ordering
        
        sseIte=sum((h(:)-h0(:)).^2);
        if sqrt(sseIte)<miniIteError 
            break 
        end
        
    end
    
    TheataRecord(:,:,t)=h;

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





% Ax+b=0 function Nest function to avoild passing many variables

    function [A,B] = PicardFdm(Value_PreviousTime)

        A=speye(length(nNode));
        B=zeros(length(nNode),1);
        
        C=theataDifFunc(h);
%         K=kFunc(h);
        K=kFieldFunc(h,Ks);

        for iz=1:nZ
            indexCenter=nodeIndex(iz);
            
            if indexCenter==0;      %is NOT a free node with index number
               continue

           %TODO this is wrong As NBC location (0 flux direction) need to be found.
            elseif    ifNbc(iz)==1  
                
                if iz==1        % a top NBC point               
                   indexUp=0;       
                   indexDown=nodeIndex(iz+1);

                   %  Forge a up ghost point 
                   nbcValue=0;  
                   hUp=  h(iz+1)- 2* deltaZ* nbcValue; %TODO NBCvalue free control
                   
                   hDown=h(iz+1);     

                   kHalfDown =(K(iz,1)+K(iz+1,1))/2;       
                   kHalfUp=kHalfDown; %approximation 
                   cCenter=C(iz,1);    
                   
                elseif iz==nZ       %if a bottom NBC point
                   indexUp=nodeIndex(iz-1);     
                   indexDown=0;              
                   
                   hUp=h(iz-1);
                   
                   % Forge a down ghost point 
                   nbcValue=-0;  
                   hDown=  h(indexUp)- 2* deltaZ* nbcValue; %TODO NBCvalue

                   kHalfUp   =(K(iz,1)+K(iz-1,1))/2;
                   kHalfDown =kHalfUp;       %approximation        
                   cCenter=C(iz,1); 
                else
                    error('innier Neumann BC point not support yet')
                end            
            else        % if Normal inner point
                                        
               indexUp=nodeIndex(iz-1);
               indexDown=nodeIndex(iz+1);
               hUp=h(iz-1);
               hDown=h(iz+1);     
               
                kHalfUp   =(K(iz,1)+K(iz-1,1))/2;
                kHalfDown =(K(iz,1)+K(iz+1,1))/2;
                cCenter=C(iz,1);
            end   

                wUp   = -kHalfUp  ./  deltaZ^2;
                wDown = -kHalfDown./  deltaZ^2;

                wCenter=cCenter/deltaTime-wUp-wDown;

                b=(kHalfDown-kHalfUp)/ deltaZ- Value_PreviousTime(iz,1)*cCenter/deltaTime;

                %modify if neighbours are DBC points           
                b=b + wUp * hUp * ~indexUp...
                    + +wDown * hDown * ~indexDown;
                
                if indexUp>0 A(indexCenter,indexUp)=wUp; end
                if indexDown>0 A(indexCenter,indexDown)=wDown; end

                A(indexCenter,indexCenter)=wCenter;
                B(indexCenter,1)=b;
        end
    end

end 



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




function Ks=permeabilityField(X)
% pointCoordinate=[X(:),Z(:)];
lenscale=10; %larger number means less stochastic field. Thus less smooth.
[nX,dimX]=size(X);

%calculate distance matrix
distance = pdist(X);
distanceMatrix = squareform(distance);

%calculate covariance matrix
covMatrix=exp(-distanceMatrix./lenscale);    

% KL decomposition on covariance matrix via SVD/eigen decomposition
% [klBasis,klEigenValue] = eigs(covMatrix,nY*nX); 
[klBasis,klEigenValue,~] = svds(covMatrix,nX); 


% [nKlBasis,~]=sizes(klBasis);


%Generate independent normal samples 
seed=100;
rng(seed);
sample= randn(nX,1);

%make multivariate Gaussian distributions with samples. zero mean.
%Covariance specified though KL basis.
Ks=klBasis*sqrt(klEigenValue)*sample;
% Ks=reshape(Ks,nY,nX);

% a log (multi) normal permeability field
Ks=exp(Ks);

end

