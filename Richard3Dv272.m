function [  ] = Richard3Dv272()
% 3D Richards equation with constant boundary condition 
% H based Richards equation
%
% First edition: Weix 25/04/2017 
%
% point type N: normal (N is natural number)
%  		  	 0: Dirichlet BC
%		  	 -N: Neuman BC 
%
% Version IV: 
%
% Version 1.40 : Weix 12/04/2017 
% improve the flexibilities (for different BC and domain) by
% introducing the indexMatrix and use 0 to indicate D-BC.
% Version 1.50 : Weix 13/04/2017 
% Vectorization
% Version 2.70: 24/04/2017 
% add Updata Neumann boundary condition. Update the way points are
% accessed. (see Richard1Dv27 for more history)
% Version 2.72: 25/04/2017 created from Richard2Dv272().
%%
tic
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
lengthTime=300;
deltaTime=1;
nTime=lengthTime/deltaTime;

% Iteration setup
nMaxIteration=1000;
miniIteError=0.1;

% Mesh
% [X,Y,Z] = meshgrid(0:deltaX:lengthX,0:deltaY:lengthY,0:deltaZ:lengthZ);

[Z,X,Y] = ndgrid(0:deltaZ:lengthZ,0:deltaX:lengthX,0:deltaY:lengthY);
% [x,y,z] = meshgrid(0:deltaX:lengthX,0:deltaY:lengthY,0:deltaZ:lengthZ);

% Permeability field

Ks=permeabilityField([Z(:),X(:),Y(:)])*0.01;
Ks=reshape(Ks,nZ,nX,nY);

% slice(Ks,Z,X,Y)

scatter3(X(:),Y(:),Z(:),Ks(:)*10000,Ks(:)*10000)

    %None of the following works
%     slice(Ks,X,Y,Z)
%     slice(Ks,Z,X,Y)
%     slice(X,Y,Z,Ks,0:deltaX:lengthX,0:deltaY:lengthY,0:deltaZ:lengthZ)
%     slice(Z,X,Y,Ks,0:deltaZ:lengthZ,0:deltaX:lengthX,0:deltaY:lengthY)

% surf(X,X,Y,Ks)

% Z2=reshape(Z(:),nZ,nX,nY);      %CHECKed correct!   

%plot permeability field
% pcolor(X,Z,Ks)
% shading interp;
% colormap jet;

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

% Option 2 Neumann BC on left and right 
    nodeIndex=zeros(nZ,nX,nY);
    nodeIndex(2:end-1,1:end,1:end)=reshape(uint32(1:(nZ-2)*nX*nY), (nZ-2),nX,nY);  
    
    nodeInFieldIndex=find(nodeIndex);

    nodeIndex(:,:,1)=-nodeIndex(:,:,1);
    nodeIndex(:,:,end)=-nodeIndex(:,:,end);
    
    nodeIndex(:,1,2:end-1)=-nodeIndex(:,1,2:end-1);
    nodeIndex(:,end,2:end-1)=-nodeIndex(:,end,2:end-1);
    
    H_init=ones(nZ,nX,nY)*-61.5;
    H_init(1,:,:)=ones(nX,nY)*-20.7;
    H_init(end,:,:)=ones(nX,nY)*-61.5;
    

%% MAIN
nNode=length(nodeIndex(nodeIndex~=0));

H=H_init;
for t=1:nTime
      
    H_PreviousTime= H;
    for k=1:nMaxIteration 
        
        H0=H;
        
        [A,B] = PicardFdm(H_PreviousTime);
        hFree=A\(-B);
                
%         H(find(nodeIndex))=hFree;       %pay extra attention to ordering
        H(nodeInFieldIndex)=hFree; 
        
        sseIte=sum((H(:)-H0(:)).^2);
        if sqrt(sseIte)<miniIteError 
            break 
        end
        
    end
    
    TheataRecord(:,:,:,t)=H;

end


toc
    
figure(1)
scatter3(X(:),Y(:),Z(:),Ks(:)*10000,Ks(:)*10000)

figure(2)
scatter3(X(:),Y(:),Z(:),abs(H(:)),abs(H(:)))

% figure(5)
% p = patch(isosurface(H,-60));
% isonormals(H,p)
% p.FaceColor = 'red';
% p.EdgeColor = 'none';
% daspect([1 1 1])
% view(3); 
% axis([nX,nY,nZ ])

figure(3)
for t=1:nTime
    clf
    Ht=abs(TheataRecord(:,:,:,t));
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
    title(sprintf('time=%i',t))
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
slice(H,nX,1:5:nZ,1)



figure(4)
for t=1:nTime
%     surf(X(:,:,sliceY),Z(:,:,sliceY),TheataRecord(:,:,sliceY,t))
%     shading interp;
%     slice(TheataRecord(:,:,:,t),nX,nY,0:nZ/4:nZ)
    slice(TheataRecord(:,:,:,t),nX-3,1:5:nZ,1)
    view(-60,40)
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;
end


figure(5)
for t=1:4:nTime
    hVector=(TheataRecord(:,:,:,t)+100)*20;    
%     scatter3(X(:),Y(:),Z(:),(hVector(:)),(Ks(:)*10000),'fill')
    scatter3(X(:),Y(:),Z(:),(hVector(:)),(hVector(:)),'fill')
%     shading interp;
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;
    
end





% Ax+b=0 function Nest function to avoild passing many variables

    function [A,B] = PicardFdm(Value_PreviousTime)

        A=speye(length(nNode));
        B=zeros(length(nNode),1);
        
        C=theataDifFunc(H);
%         K=kFunc(H);
        K=kFieldFunc(H,Ks);

        for iZ=1:nZ
            for iX=1:nX
                for iY=1:nY
            
                    indexCenter=nodeIndex(iZ,iX,iY);

                    switch sign(indexCenter)                    
                        case 0      %is NOT a free node with index number
                            continue
                        case -1     %TODO this is wrong As NBC location (0 flux direction) need to be found.

                            %---------------------------------------------
                            if iZ==1        % a top NBC point               
                               indexUp=0;                             
                               %Forge a up ghost point 
                               nbcValue=0;  
                               hUp=  H(iZ+1,iX,iY)- 2* deltaZ* nbcValue; %TODO NBCvalue free control

                               kHalfDown =(K(iZ,iX,iY)+K(iZ+1,iX,iY))/2;                         
                               kHalfUp=kHalfDown; %Vitual K  
                            else
                               indexUp=nodeIndex(iZ-1,iX,iY);
                               kHalfUp  =(K(iZ,iX,iY)+K(iZ-1,iX,iY))/2; 
                               hUp=H(iZ-1,iX,iY);
                            end

                            if iZ==nZ       %if a bottom NBC point             
                               indexDown=0;                             
                               %Forge a up ghost point 
                               nbcValue=0;  
                               hDown=  H(iZ-1,iX,iY)- 2* deltaZ* nbcValue; %TODO NBCvalue free control

                               kHalfUp   =(K(iZ,iX,iY)+K(iZ-1,iX,iY))/2;
                               kHalfDown =kHalfUp;       %Vitual K     
                            else 
                               indexDown=nodeIndex(iZ+1,iX,iY);
                               kHalfDown  =(K(iZ,iX,iY)+K(iZ+1,iX,iY))/2; 
                               hDown=H(iZ+1,iX,iY);  
                            end
                            
                            %---------------------------------------------
                            if iX==1        % a left NBC point
                               indexLeft=0;                             
                               %Forge a up ghost point 
                               nbcValue=0;  
                               hLeft=  H(iZ,iX+1,iY)- 2* deltaX* nbcValue; %TODO NBCvalue free control

                               kHalfRight =(K(iZ,iX,iY)+K(iZ,iX+1,iY))/2;                         
                               kHalfLeft=kHalfRight; %Vitual K  

                            else
                               indexLeft=nodeIndex(iZ,iX-1,iY);
                               kHalfLeft  =(K(iZ,iX,iY)+K(iZ,iX-1,iY))/2; 
                               hLeft=H(iZ,iX-1,iY);
                            end

                            if iX==nX        % a right NBC point               
                               indexRight=0;                             
                               %Forge a up ghost point 
                               nbcValue=0;  
                               hRight=  H(iZ,iX-1,iY)- 2* deltaX* nbcValue; %TODO NBCvalue free control

                               kHalfLeft =(K(iZ,iX,iY)+K(iZ,iX-1,iY))/2;                         
                               kHalfRight=kHalfLeft; %Vitual K  

                            else
                               indexRight=nodeIndex(iZ,iX+1,iY);
                               kHalfRight  =(K(iZ,iX,iY)+K(iZ,iX+1,iY))/2; 
                               hRight=H(iZ,iX+1,iY);
                            end
                            
                            %---------------------------------------------
                            if iY==1        % a front NBC point
                               indexFront=0;                             
                               %Forge a up ghost point 
                               nbcValue=0;  
                               hFront=  H(iZ,iX,iY+1)- 2* deltaY* nbcValue; %TODO NBCvalue free control

                               kHalfBack =(K(iZ,iX,iY)+K(iZ,iX,iY+1))/2;                         
                               kHalfFront=kHalfBack; %Vitual K  

                            else
                               indexFront=nodeIndex(iZ,iX,iY-1);
                               kHalfFront  =(K(iZ,iX,iY)+K(iZ,iX,iY-1))/2; 
                               hFront=H(iZ,iX,iY-1);
                            end
                            
                            if iY==nY        % a front NBC point
                               indexBack=0;                             
                               %Forge a up ghost point 
                               nbcValue=0;  
                               hBack=  H(iZ,iX,iY-1)- 2* deltaY* nbcValue; %TODO NBCvalue free control

                               kHalfFront =(K(iZ,iX,iY)+K(iZ,iX,iY-1))/2;                         
                               kHalfBack=kHalfFront; %Vitual K  

                            else
                               indexBack=nodeIndex(iZ,iX,iY+1);
                               kHalfBack =(K(iZ,iX,iY)+K(iZ,iX,iY+1))/2; 
                               hBack=H(iZ,iX,iY+1);
                            end
                            %---------------------------------------------

                        case 1      % if Normal inner point

                            indexUp=nodeIndex(iZ-1,iX,iY);
                            indexDown=nodeIndex(iZ+1,iX,iY);
                            indexLeft=nodeIndex(iZ,iX-1,iY);
                            indexRight=nodeIndex(iZ,iX+1,iY);                           
                            indexFront=nodeIndex(iZ,iX,iY-1);
                            indexBack=nodeIndex(iZ,iX,iY+1);
                            
                            hUp=H(iZ-1,iX,iY);
                            hDown=H(iZ+1,iX,iY);
                            hLeft=H(iZ,iX-1,iY);
                            hRight=H(iZ,iX+1,iY);
                            hFront=H(iZ,iX,iY-1);
                            hBack=H(iZ,iX,iY+1);
                            
                            kHalfUp   =(K(iZ,iX,iY)+K(iZ-1,iX,iY))/2;
                            kHalfDown =(K(iZ,iX,iY)+K(iZ+1,iX,iY))/2;
                            kHalfLeft =(K(iZ,iX,iY)+K(iZ,iX-1,iY))/2;
                            kHalfRight=(K(iZ,iX,iY)+K(iZ,iX+1,iY))/2;                                                        
                            kHalfFront=(K(iZ,iX,iY)+K(iZ,iX,iY-1))/2;   
                            kHalfBack =(K(iZ,iX,iY)+K(iZ,iX,iY+1))/2;  
                            
                        otherwise 
                            error('unknown node type');

                    end   

                    cCenter=C(iZ,1);   

                    wUp   = -kHalfUp  ./deltaZ^2;
                    wDown = -kHalfDown./deltaZ^2;
                    wLeft = -kHalfLeft./deltaX^2;
                    wRight= -kHalfRight./deltaX^2;                    
                    wFront= -kHalfFront./deltaY^2;
                    wBack = -kHalfBack ./deltaY^2;
                    
                    
                    

                    wCenter=cCenter/deltaTime-wUp-wDown-wLeft-wRight-wFront-wBack;

                    b=(kHalfDown-kHalfUp)/deltaZ-H_PreviousTime(iZ,iX,iY)*cCenter/deltaTime;

                    %modify if neighbours are DBC points           
                    b=b + wUp    * hUp    * ~indexUp...
                        + wDown  * hDown  * ~indexDown...
                        + wLeft  * hLeft  * ~indexLeft...
                        + wRight * hRight * ~indexRight...
                        + wFront * hFront * ~indexFront...
                        + wBack  * hBack  * ~indexBack;


                    indexUp=abs(indexUp);
                    indexDown=abs(indexDown);
                    indexLeft=abs(indexLeft);
                    indexRight=abs(indexRight);
                    indexFront=abs(indexFront);
                    indexBack=abs(indexBack);
                    
                    indexCenter=abs(indexCenter);

                    if indexUp>0 A(indexCenter,indexUp)=wUp; end
                    if indexDown>0 A(indexCenter,indexDown)=wDown; end
                    if indexLeft>0 A(indexCenter,indexLeft)=wLeft; end
                    if indexRight>0 A(indexCenter,indexRight)=wRight; end
                    if indexFront>0 A(indexCenter,indexFront)=wFront; end
                    if indexBack>0 A(indexCenter,indexBack)=wBack; end                    
                               
                    A(indexCenter,indexCenter)=wCenter;
                    B(indexCenter,1)=b;
                end %iY loop
            end %iX loop
        end %iZ loop
    end

end 



function theata=theataFunc(H)
theataS=0.287;
theataR=0.075;
alpha=1.611e6;
beta=3.96;

result=alpha.*(theataS-theataR)/(alpha+abs(H).^beta)+theataR;
end

function theataDif=theataDifFunc(H)
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

theataDif=-alpha.*(theata_s-theata_r).*-1.*(alpha+abs(H).^beta).^(-2).*abs(H).^(beta-1);

end

function result=kFunc(H)
rho=1.175e6;
r=4.74;
k_s=0.00944;

result=k_s.*rho./(rho+abs(H).^r);
end

function result=kFieldFunc(H,ks)
% H and k must be the same sizes
rho=1.175e6;
r=4.74;

result=ks.*rho./(rho+abs(H).^r);
end




function Ks=permeabilityField(X)
% pointCoordinate=[X(:),Z(:)];
lenscale=100; %larger number means less stochastic field. Thus less smooth.
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

