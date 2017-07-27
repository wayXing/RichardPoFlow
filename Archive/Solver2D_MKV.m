function [  ] = Solver2D_MKV()
% 2D Richards equation with constant boundary condition 
% h based Richards equation
%
% Version IV: 
%
% Version 1.40 : Weix 12/04/2017 
% improve the flexibilities (for different BC and domain) by
% introducing the indexMatrix and use 0 to indicate D-BC.
% Version 1.50 : Weix 13/04/2017 
% Vectorization
%%
tic
% nNode=10;
% nTime=100;

% Spatial setup
lengthZ=40;
deltaZ=1;
nNodeZ=lengthZ/deltaZ-1;

lengthX=40;
deltaX=2;
nNodeX=lengthX/deltaX-1;

% Temporal setup
lengthTime=300;
deltaTime=1;
nTime=lengthTime/deltaTime;

% Iteration setup
nMaxIteration=1000;
miniIteError=0.1;


% Mesh
[X,Z] = meshgrid(0:deltaX:lengthX,0:deltaZ:lengthZ);

Ks=permeabilityField(X,Z)*0.01;
pcolor(Ks)
shading interp;
colormap jet;

% integer shows unknowns and 0 shows DBC
nUnknown=nNodeZ*nNodeX;
indexMatrix=zeros(nNodeZ+2, nNodeX+2);
indexMatrix(2:end-1,2:end-1)=reshape(uint32(1:nUnknown), (nNodeZ), (nNodeX));

% indexUnknown = reshape(uint32(1:nUnknown), (nNodeZ-2), (nNodeX-2));

% indexMatrix=zeros(nNodeZ, nNodeX);
% indexMatrix(2:end-1,2:end-1)=indexUnknown;



H=zeros(nNodeZ,nNodeX);
% H_temp=zeros(nNode,nMaxIteration);

%%% 
%initial state
H_init=zeros(nNodeZ,nNodeX);
H_init(:)=-61.5;

%BC
% bcLeft=ones(nNodeZ,1)*-20.7;
% bcRight=ones(nNodeZ,1)*-61.5;
% bcTop=ones(nNodeX,1)*-20.7;
% bcBottom=ones(nNodeX,1)*-61.5;


    % A more interesting setup
    bcLeft=ones(nNodeZ,1)*-20.7;
    bcRight=ones(nNodeZ,1)*-20.7;
    bcTop=ones(nNodeX,1)*-20.7;
    bcBottom=ones(nNodeX,1)*-24.7;

    
    
H_all=zeros(nNodeZ+2,nNodeX+2);
H_all(1,2:end-1)=bcTop;
H_all(end,2:end-1)=bcBottom;

H_all(2:end-1,1)=bcLeft;
H_all(2:end-1,end)=bcRight;
H_all(2:end-1,2:end-1)=H_init;


% update BC to initial field
% H_init(:,1)=bcLeft;
% H_init(:,end)=bcRight;
% 
% H_init(end,:)=bcBottom;
% H_init(1,:)=bcTop;

%%
% indicate the BC points
temp=zeros(nNodeZ,nNodeX);
% temp(1,:)=1;
% ifTopDbc=temp(:);   %if top neighbour a constant B.C. 
temp(1,:)=bcTop;      %this structure might be convient but require storage
topDbcValue=temp(:);

temp=zeros(nNodeZ,nNodeX);
% temp(end,:)=1;
% ifBottomDbc=temp(:);   
temp(end,:)=bcBottom;
bottomDbcValue=temp(:);

temp=zeros(nNodeZ,nNodeX);
% temp(:,1)=1;
% ifLeftDbc=temp(:);  
temp(:,1)=bcLeft;
leftDbcValue=temp(:);

temp=zeros(nNodeZ,nNodeX);
% temp(:,end)=1;
% ifRightDbc=temp(:);   
temp(:,end)=bcRight;
rightDbcValue=temp(:);

%% Another way. We can just use self organize as this need only onces.
% TopDbcValue=[zeros(1, size(H_all,2)),H_all];




%% MAIN
TheataRecord(:,:,1)=H_init;
H=H_init;

% C=ones(nNodeZ,nNodeX)*1234567;

for t=1:nTime
    
    H_PreviousTime=H;
    
    
    for k=1:nMaxIteration 
    
        H0=H;
        
%         C=theataDifFunc(H);
%         K=kFunc(H);
%         K=kFieldFunc(H,Ks);
        
%         D=D_Func(H);
        
        %initialize
        A=zeros(nUnknown);
        B=zeros(nUnknown,1);
        

        %% New Vectorization  %Assemble Ax+B=0;
%         [I,J] = ind2sub([nNodeZ,nNodeX],1:nUnknown)
        
%         C=theataDifFunc(H);
%         C_bcTop=theataDifFunc(bcTop);
%         C_bcBottom=theataDifFunc(bcBottom);
%         C_bcLeft=theataDifFunc(bcLeft);
%         C_bcRight=theataDifFunc(bcRight);
%         
%         
%         K=kFunc(H);
%         K_bcTop=kFunc(bcTop);
%         K_bcBottom=kFunc(bcBottom);
%         K_bcLeft=kFunc(bcLeft);
%         K_bcRight=kFunc(bcRight);
        
%         cField=theataDifFunc(H);
%         kField=kFunc(H);
        
        C_all=theataDifFunc(H_all);
        K_all=kFunc(H_all);
        K_all=kFieldFunc(H_all,Ks);
        
        
        C=C_all(2:end-1,2:end-1);
        K=K_all(2:end-1,2:end-1);
        
        
%         zdiffC_all=diff(C_all,1,1);
%         xdiffC_all=diff(C_all,1,2);
        
        zdiffK_all=diff(K_all,1,1);
        xdiffK_all=diff(K_all,1,2);
        
        
        wUp   = -1./deltaZ^2 .* (K - zdiffK_all(1:end-1,2:end-1)./2);
        wDown = -1./deltaZ^2 .* (K + zdiffK_all(2:end,  2:end-1)./2);
        
        wLeft = -1./deltaX^2 .* (K - xdiffK_all(2:end-1,1:end-1)./2);
        wRight= -1./deltaX^2 .* (K + xdiffK_all(2:end-1,2:end)./2);
        
        wCenter=C./deltaTime-wUp-wDown-wLeft-wRight;
        
        b= (zdiffK_all(2:end, 2:end-1)+ zdiffK_all(1:end-1,2:end-1))./2 ./deltaZ ...
           -H_PreviousTime .* C ./ deltaTime;
        
        
        %update BC neighbour 
      
%         b2= b(:)+ wUp(:)   .* ~(topDbcValue==0).*topDbcValue...
%                 + wDown(:) .* ~(bottomDbcValue==0).*bottomDbcValue...
%                 + wLeft(:) .* ~(leftDbcValue==0).*leftDbcValue...
%                 + wRight(:).* ~(rightDbcValue==0).*rightDbcValue;
            
        B= b(:)+ wUp(:)   .*topDbcValue...
                + wDown(:) .*bottomDbcValue...
                + wLeft(:) .*leftDbcValue...
                + wRight(:).*rightDbcValue;    
            
         
        wUp   = wUp(:)   .* (topDbcValue==0); 
        wDown = wDown(:) .* (bottomDbcValue==0); 
        wLeft = wLeft(:) .* (leftDbcValue==0); 
        wRight= wRight(:).* (rightDbcValue==0); 
        
        wCenter=wCenter(:);
        
        A=  diag(wUp(2:end),1)+diag(wDown(1:end-1),-1)...
            +diag(wLeft(nNodeZ+1:end),-nNodeZ)+diag(wRight(1:end-nNodeZ),nNodeZ)...       
            +diag(wCenter,0);

        
%         %trime for BC
%         indexIfKnown=~any(A);
%         indexIfUnknown=any(A);
% 
%         A=A(indexIfUnknown,indexIfUnknown);
%         B=B(indexIfUnknown);

        %% 
        h=A\(-B);         
        H=reshape(h,nNodeZ,nNodeX);
        
        
        sseIte=sum((H(:)-H0(:)).^2);
        if sqrt(sseIte)<miniIteError 
            break 
        end
        
        H_all(2:end-1,2:end-1)=H;
        
        
    end
    
    
    
    
    TheataRecord(:,:,t)=H;
    
end
    
toc
    
figure(1)
surf(H_init);

for t=1:nTime
    surf(TheataRecord(:,:,t))
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;
    
end
    

end


function []=axbFunc(H,indexMatrix)

    [nNodeZ,nNodeX]=size(indexMatrix);
    
    C=theataDifFunc(H);
    K=kFunc(H);
    
    nUnknown=sum(any(indexMatrix));
    
    %initialize
    A=zeros(nUnknown);
    B=zeros(nUnknown,1);

    [I,J] = ind2sub(siz,IND)
    
    for i=1:nUnknown
        
        
        
    end
    
end








function [k,j]=index2kj(index)
    
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

function result=kFieldFunc(h,k)
% h and k must be the same sizes
rho=1.175e6;
r=4.74;

result=k.*rho./(rho+abs(h).^r);
end



function Ks=permeabilityField(X,Z)
% pointCoordinate=[X(:),Z(:)];
lenscale=10; %larger number means less stochastic field. Thus less smooth.

[nY,nX]=size(X);

distance = pdist([X(:),Z(:)]);
distanceMatrix = squareform(distance);

covMatrix=exp(-distanceMatrix./lenscale);    

% [klBasis,klEigenValue] = eigs(covMatrix,nY*nX); 
[klBasis,klEigenValue,~] = svds(covMatrix,nY*nX); 


% [nKlBasis,~]=sizes(klBasis);

% seed=100;
% rng(seed);
sample= randn(nY*nX,1);


Ks=klBasis*sqrt(klEigenValue)*sample;
Ks=reshape(Ks,nY,nX);

Ks=exp(Ks);

end