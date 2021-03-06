function [Ar,Br]=picardAxbFormDr(mesh,previousZh,deltaT,Vh,Zh,Vk,Zk,Vc,Zc)
%calculate linear system equation Ax=B using Picards scheme with dimension
%reduction DEMI.
%This is not a practical function but very useful for understanding the
%program.
%
% Input parameters:
%   mesh             -mseh structure 
%   previousH        -value of mesh.H at last time step
% Output parameters:
%   A,B              -A*mesh.H_new=B;
%
% Examples: see Demo
%
% % Short description of example, followed by Matlab code line
% >> matlab code of example
%
% See also: 
% Author:   Wei Xing
% History:  10/05/2017  file created
%
%
%     Kr=Vk*Zk; %approximation to K
%     Cr=Vc*Zc;
% 
%     C=mesh.C;
%     K=mesh.K;


%%  Auxiliary variable   
dbcIndex=find(mesh.dbcFlag);     %specify DBC index for later fitting in value
nodeIndex=find(~mesh.dbcFlag);   %specify free node index
nNode=sum(~mesh.dbcFlag);        %number of free node   

deltaZ=mesh.deltaZ;
nZ=mesh.nZ;

% P=diag(mesh.dbcFlag);              %Picking up matrix. maybe useful
% P=spdiags(mesh.dbcFlag,0,nZ,nZ)    %Picking up matrix. maybe useful

% P=diag(mesh.dbcFlag);     

Pwave=spdiags(mesh.dbcFlag,0,nZ,nZ);    %Picking up matrix. maybe useful
P    =spdiags(~mesh.dbcFlag,0,nZ,nZ);    %Picking up matrix. maybe useful

Pwave( ~any(Pwave,2), : ) = [];  %remove zero rows
P(~any(P,2), : )          = [];  %remove zero rows

% PP=P'*P;
% PPwave=Pwave'*Pwave;
% 
% VhWave=Pwave'*Pwave*Vh;
% Vh=P'*P*Vh;


%% form shift matrix. This matrix is not always needed, Depends on the method.
iMethod=1;
switch iMethod
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
    
% KrLowShift=lowShift1Eye*Vk; 
% KrUpShift=UpShift1Eye*Vk; 
% Kr=Vk*Zk;     %approximation to K using Kr


%% 
%may be useful
% for i=1:length(Zk)
%     centerDiagR(:,i) = (2.*Vk(:,i)+ lowShift1Eye*Vk(:,i)+UpShift1Eye*Vk(:,i)) *Zk(i)./(2*deltaZ^2)+ Vc(:,i)*Zc(i)./deltaT;
%     upDiagR(:,i)     = (Vk(:,i) + lowShift1Eye*Vk(:,i)).*Zk(i)./(-2*deltaZ^2);  
%     downDiagR(:,i)   = (Vk(:,i) + UpShift1Eye*Vk(:,i)) .*Zk(i)./(-2*deltaZ^2);  
% end
% 
%     centerDiagR = (2.*Vk+ lowShift1Eye*Vk+UpShift1Eye*Vk(:,i)) ./(2*deltaZ^2).*Zk+ Vc./deltaT *Zc;
%     upDiagR     = (Vk + lowShift1Eye*Vk)./(-2*deltaZ^2) *Zk;  
%     downDiagR   = (Vk + UpShift1Eye*Vk) ./(-2*deltaZ^2) *Zk;  

    
%     Ar=[];
    [~,nDimR]=size(Vh);
    Ar=zeros(nDimR,nDimR);
    Ar2=zeros(nDimR,nDimR);
    
    for i=1:length(Zk)
        %loop can be replace using tensor product.see note.
       
%         ArC=Ar +Vh'* spdiags(centerDiagR(:,i),0,nZ,nZ)*Vh .*Zk(i);  
%                +Vh'* spdiags(upDiagR(:,i),    0,nZ,nZ)*lowShift1Eye *Vh 
%                +Vh'* spdiags(downDiagR(:,i),  0,nZ,nZ)*UpShift1Eye  *Vh;  

        %ArV means Ar's basis V. a matrix basis.
        centerDiagRVPart1(:,i)= (2.*Vk(:,i)+ lowShift1Eye*Vk(:,i)+UpShift1Eye*Vk(:,i))./(2*deltaZ^2);
        centerDiagRVPart2(:,i)= Vc(:,i)./deltaT;
        upDiagR(:,i)          = (Vk(:,i)   + lowShift1Eye*Vk(:,i))                    ./(-2*deltaZ^2);
        downDiagR(:,i)        = (Vk(:,i)   + UpShift1Eye*Vk(:,i))                     ./(-2*deltaZ^2);
        
%         ArV(:,:,i)=  Vh'*spdiags(centerDiagRVPart1(:,i),0,nZ,nZ)              *Vh .*Zk(i)...
%                     +Vh'*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)              *Vh .*Zc(i)...    
%                     +Vh'*spdiags(upDiagR(:,i)          ,0,nZ,nZ)*lowShift1Eye *Vh .*Zk(i)...
%                     +Vh'*spdiags(downDiagR(:,i)        ,0,nZ,nZ)*UpShift1Eye  *Vh .*Zk(i);
        
        %Add considerations to BC        
        ArV(:,:,i)=  Vh'*P'*P*spdiags(centerDiagRVPart1(:,i),0,nZ,nZ)              *P'*P*Vh .*Zk(i)...
                    +Vh'*P'*P*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)              *P'*P*Vh .*Zc(i)...    
                    +Vh'*P'*P*spdiags(upDiagR(:,i)          ,0,nZ,nZ)*lowShift1Eye *P'*P*Vh .*Zk(i)...
                    +Vh'*P'*P*spdiags(downDiagR(:,i)        ,0,nZ,nZ)*UpShift1Eye  *P'*P*Vh .*Zk(i);                
                           
        Ar=Ar+ArV(:,:,i);
        
        %%TODO consider very carefully whether to use P(pickup matrix) to take
        %%known node.
        %if so, use follow
%         ArV(i)=  (P*Vh)'*(P*spdiags(centerDiagRVPart1(i),0,nZ,nZ)             *P') *(P*Vh) .*Zk(i)...
%                 +(P*Vh)'*(P*spdiags(centerDiagRVPart2(i),0,nZ,nZ)             *P') *(P*Vh) .*Zc(i)...    
%                 +(P*Vh)'*(P*spdiags(upDiagR(i)          ,0,nZ,nZ)*lowShift1Eye*P') *(P*Vh) .*Zk(i)...
%                 +(P*Vh)'*(P*spdiags(downDiagR(i)        ,0,nZ,nZ)*UpShift1Eye *P') *(P*Vh) .*Zk(i);    

        ArV2(:,:,i)= Vh'*P'*P*spdiags(centerDiagRVPart1(:,i),0,nZ,nZ)              *Pwave'*Pwave*Vh .*Zk(i)...
                    +Vh'*P'*P*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)              *Pwave'*Pwave*Vh .*Zc(i)...    
                    +Vh'*P'*P*spdiags(upDiagR(:,i)          ,0,nZ,nZ)*lowShift1Eye *Pwave'*Pwave*Vh .*Zk(i)...
                    +Vh'*P'*P*spdiags(downDiagR(:,i)        ,0,nZ,nZ)*UpShift1Eye  *Pwave'*Pwave*Vh .*Zk(i);                
                         
        Ar2=Ar2+ArV2(:,:,i);
        
    end
    
    
    
    %% form B vector 
    
    iMethod=1;
    switch iMethod 
        case 1
        Br =   Vh'*P'*P * -(UpShift1Eye-lowShift1Eye)./(2*deltaZ)*Vk*Zk...
              +Vh'*P'*P *   ((Vh*previousZh) ./deltaT .*(Vc *Zc));
      
        %case 2 representation is needed for real ROM
        case 2
        Brk=Vh'*P'*P * -(UpShift1Eye-lowShift1Eye)./(2*deltaZ)*Vk*Zk;
        Br=Brk;
        for i=1:size(Vh,2)
           Brhc(:,:,i)=Vh'*P'*P*spdiags(Vh(:,i),0,nZ,nZ)*Vc *previousZh(i)*Zc./deltaT;
           Br         = Br +Brhc(:,:,i);
        end
    end
    
      
      
    %%Take away BC point in ROM  
    Br =   Br- Ar2*Zh;
    
    
% centerDiagR =@(Zk) (2.*Vk*Zk+ KrLowShift*Zk+KrUpShift*Zk)./(2*deltaZ^2)+Vc*Zc./deltaT; %R means approximation from reduced basis
% upDiagR     =@(Zk) (Vk*Zk+ KrLowShift*Zk)/(-2*deltaZ^2);  
% downDiagR   =@(Zk) (Vk*Zk+ KrUpShift*Zk)/(-2*deltaZ^2);   


% centerDiag = (2.*Vk*Zk+ KrLowShift*Zk+KrUpShift*Zk)./(2*deltaZ^2)+Vc*Zc./deltaT;       %A_center diagonal %first and last elements are meaningless                    
% upDiag     = (Vk*Zk+ KrLowShift*Zk)/(-2*deltaZ^2);                                %A_up     diagonal %first and last elements are meaningless                  
% downDiag   = (Vk*Zk+ KrUpShift*Zk)/(-2*deltaZ^2);   

%      B=B(nodeIndex)-A_all(nodeIndex,dbcIndex)*mesh.H(dbcIndex);
    %                     A=A_all(nodeIndex,nodeIndex);


%             %% form the sparse band   
%             iMethod=1;
%             switch iMethod
%                 case 1  
%                     %write 3 diagonal band  %very fast. as avoid matrix 
%                     %calculation and use a element wise operation.
%                     centerDiag = (2.*K+ lowShift1Eye*K+UpShift1Eye*K)/(2*deltaZ^2)+C/deltaT;       %A_center diagonal %first and last elements are meaningless                    
%                     upDiag     = (K+ lowShift1Eye*K)/(-2*deltaZ^2);                                %A_up     diagonal %first and last elements are meaningless                  
%                     downDiag   = (K+ UpShift1Eye*K)/(-2*deltaZ^2);                                 %A_down   diagonal %first and last elements are meaningless    
% 
%                     centerDiag = (2.*Vk*Zk+ KrLowShift*Zk+KrUpShift*Zk)./(2*deltaZ^2)+Vc*Zc./deltaT;       %A_center diagonal %first and last elements are meaningless                    
%                     upDiag     = (Vk*Zk+ KrLowShift*Zk)/(-2*deltaZ^2);                                %A_up     diagonal %first and last elements are meaningless                  
%                     downDiag   = (Vk*Zk+ KrUpShift*Zk)/(-2*deltaZ^2);   
% 
%                 case 2     
%                     %rewrite 3 diagonal band  %fast. do matrix calculation though
%                     %using sparse operations.
%                     centerDiag = (2.*speye(nZ)+lowShift1Eye+UpShift1Eye)./(2*deltaZ^2)*K+C/deltaT;   
%                     upDiag     = (speye(nZ)+lowShift1Eye)./(2*deltaZ^2)*K;
%                     downDiag   = (speye(nZ)+UpShift1Eye)./(2*deltaZ^2)*K;
%                 case 3
%                     centerDiag = (2.*speye(nZ)+circshift(speye(nZ),[1,0])+circshift(speye(nZ),[-1,0]))./(2*deltaZ^2)*K+C/deltaT;   
%                     upDiag     = (speye(nZ)+circshift(speye(nZ),[1,0]))./(2*deltaZ^2)*K;
%                     downDiag   = (speye(nZ)+circshift(speye(nZ),[-1,0]))./(2*deltaZ^2)*K;
% 
%             end
% 
%             %% make spare metrix A from bands                  
%             %Todo this part may be improved using better sparse diag
%             iMethod=4;
%             switch iMethod 
%                 case 1      %Very fast
%                     A_all=spdiags(centerDiag,0,nZ,nZ) +circshift(spdiags(upDiag,0,nZ,nZ),[0,-1]) +circshift(spdiags(downDiag,0,nZ,nZ),[0,1]);                                          
%                 case 2      %Slow. Always NOT creat full matrix and then sparse it.
%                     A_all=sparse (diag(centerDiag,0) +circshift(diag(upDiag,0),[0,-1]) +circshift(diag(downDiag,0),[0,1]));   
%                 case 3      %not necessary but show a more clear formulation
%                     left1Eye=circshift(spdiags(ones(nZ),0,nZ,nZ),[0,-1]);
%                     right1Eye=circshift(spdiags(ones(nZ),0,nZ,nZ),[0,1]);
%                     A_all=spdiags(centerDiag,0,nZ,nZ) +spdiags(upDiag,0,nZ,nZ)*left1Eye + spdiags(downDiag,0,nZ,nZ)*right1Eye;   
%                 case 4
%                     A_all=spdiags(centerDiag,0,nZ,nZ) +spdiags(upDiag,0,nZ,nZ)*lowShift1Eye + spdiags(downDiag,0,nZ,nZ)*UpShift1Eye;  
% 
% 
%                 case 5
%             %             A_all=spdiags(centerDiag,0,nZ,nZ) +spdiags(upDiag,0,nZ,nZ)*lowShift1Eye + spdiags(downDiag,0,nZ,nZ)*UpShift1Eye;  
%                 case 6
%                     A_all=centerDiag*speye(nZ) +upDiag*speye(nZ)*lowShift1Eye + downDiag*speye(nZ)*UpShift1Eye;  
%                 otherwise 
%             end
% 
%             %% form B vector 
%             B          = -(UpShift1Eye*K-lowShift1Eye*K)/(2*deltaZ)+previousH.*C/deltaT;
% 
%                 %pick free node and componsate for dbc involved 
% 
%             %% Picking up the unknown free node
%             iMethod=1;
%             switch iMethod
%                 case 1    
%                     B=B(nodeIndex)-A_all(nodeIndex,dbcIndex)*mesh.H(dbcIndex);
%                     A=A_all(nodeIndex,nodeIndex);
%                 case 2
%                     P=spdiags(dbcFlag,0,nZ,nZ);    %Picking up matrix. maybe useful
% 
%                     B=P*B-P*A_all*P'* P*mesh.H;
%                     A=P*A_all*P';
%             end
    
    
    
    
end


