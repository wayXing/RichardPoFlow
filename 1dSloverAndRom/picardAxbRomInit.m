function [romMesh]=picardAxbRomInit(mesh,Vh,VdK,Pk,VdC,Pc)
%Initilize Picard iteration using DEIM ROM
%on the 1d h-based Richards equation
%
% Input parameters:
%   mesh             -mseh structure 
%   Vh               -basis for presure head.
%   VdK              -deim basis for k term using Deim. normally called Dk and 
%                     calculted by Dk=Vk*inv(Pk'*Vk); where Vk is the pod
%                     basis and Pk is the pick up matrix 
%   Pk               -Pickup matrix comes with VdK using Deim
%   VdC              -deim basis for c term using Deim. normally called Dc and 
%                     calculted by Dc=Vc*inv(Pc'*Vc); where Vc is the pod
%                     basis and Pc is the pick up matrix 
%   Pc               -Pickup matrix comes with VdK using Deim
% Output parameters:
%  romMesh           -reduced order mesh with its special structure
%
% Examples: see Demo
%
% See also: 
% Author:   Wei Xing
% History:  22/05/2017  file created


%%  Auxiliary variable   
dbcIndex=find(mesh.dbcFlag);     %specify DBC index for later fitting in value
nodeIndex=find(~mesh.dbcFlag);   %specify free node index
nNode=sum(~mesh.dbcFlag);        %number of free node   

deltaZ=mesh.deltaZ;
nZ=mesh.nZ;

Pwave=spdiags(mesh.dbcFlag,0,nZ,nZ);     %Picking up matrix of boundary condition
P    =spdiags(~mesh.dbcFlag,0,nZ,nZ);    %Picking up matrix of free nodel

Pwave( ~any(Pwave,2), : ) = [];  %remove zero rows
P(~any(P,2), : )          = [];  %remove zero rows

nDeimK=size(VdK,2);
nDeimC=size(VdC,2);
nPod  =size(Vh,2);

Pk=sparse(logical(Pk));
Pc=sparse(logical(Pc));
%%


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
    
    
%     Ar=zeros(nPod,nPod);
%     ArBC=zeros(nPod,nPod);
    
    %k related terms
    for i=1:nDeimK  %ONLY allow for nDeimK=nDeimC
        %loop can be replaced using tensor product.see note.
       
        %ArV means Ar's basis V. a matrix basis.
        centerDiagRVPart1(:,i)= (2.*VdK(:,i)+ lowShift1Eye*VdK(:,i)+UpShift1Eye*VdK(:,i))./(2*deltaZ^2);
%         centerDiagRVPart2(:,i)= Vc(:,i);
        upDiagR(:,i)          = (VdK(:,i)   + lowShift1Eye*VdK(:,i))                    ./(-2*deltaZ^2);
        downDiagR(:,i)        = (VdK(:,i)   + UpShift1Eye*VdK(:,i))                     ./(-2*deltaZ^2);
          
        %Ar matrix part 1. It consists of (length(Zk)) basis matrix. 
%         ArV(:,:,i)=  Vh'*P'*P*spdiags(centerDiagRVPart1(:,i),0,nZ,nZ)              *P'*P*Vh .*Zk(i)...
%                     +Vh'*P'*P*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)              *P'*P*Vh .*Zc(i)...    
%                     +Vh'*P'*P*spdiags(upDiagR(:,i)          ,0,nZ,nZ)*lowShift1Eye *P'*P*Vh .*Zk(i)...
%                     +Vh'*P'*P*spdiags(downDiagR(:,i)        ,0,nZ,nZ)*UpShift1Eye  *P'*P*Vh .*Zk(i);  
                
        Ark(:,:,i) = Vh'*P'*P*spdiags(centerDiagRVPart1(:,i),0,nZ,nZ)              *P'*P*Vh...       
                    +Vh'*P'*P*spdiags(upDiagR(:,i)          ,0,nZ,nZ)*lowShift1Eye *P'*P*Vh...
                    +Vh'*P'*P*spdiags(downDiagR(:,i)        ,0,nZ,nZ)*UpShift1Eye  *P'*P*Vh;    %Ar k related part
%         Arc(:,:,i) = Vh'*P'*P*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)              *P'*P*Vh;    %Ar c related part     
                       
        %Ar=Ar+Ark(:,:,i).*Zk(i) +Arc(:,:,i).*Zc(i)./deltaT;
                
        %consider very carefully How to deal with BC point in the ROM
%         ArV2(:,:,i)= Vh'*P'*P*spdiags(centerDiagRVPart1(:,i),0,nZ,nZ)              *Pwave'*Pwave*Vh .*Zk(i)...
%                     +Vh'*P'*P*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)              *Pwave'*Pwave*Vh .*Zc(i)...    
%                     +Vh'*P'*P*spdiags(upDiagR(:,i)          ,0,nZ,nZ)*lowShift1Eye *Pwave'*Pwave*Vh .*Zk(i)...
%                     +Vh'*P'*P*spdiags(downDiagR(:,i)        ,0,nZ,nZ)*UpShift1Eye  *Pwave'*Pwave*Vh .*Zk(i);                
        
        ArBCk(:,:,i) = Vh'*P'*P*spdiags(centerDiagRVPart1(:,i),0,nZ,nZ)              *Pwave'*Pwave*Vh...       
                      +Vh'*P'*P*spdiags(upDiagR(:,i)          ,0,nZ,nZ)*lowShift1Eye *Pwave'*Pwave*Vh...
                      +Vh'*P'*P*spdiags(downDiagR(:,i)        ,0,nZ,nZ)*UpShift1Eye  *Pwave'*Pwave*Vh;      %ArBC k related part
%         ArBCc(:,:,i) = Vh'*P'*P*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)              *Pwave'*Pwave*Vh;      %ArBC c related part 
                
        %ArBC=Ar+ArBCk(:,:,i).*Zk(i) +ArBCc(:,:,i).*Zc(i)        
        
    end
    
    %c related terms
    for i=1:nDeimC
        centerDiagRVPart2(:,i)= VdC(:,i);
        Arc(:,:,i) = Vh'*P'*P*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)              *P'*P*Vh;              %Ar   c related part     
        ArBCc(:,:,i) = Vh'*P'*P*spdiags(centerDiagRVPart2(:,i),0,nZ,nZ)            *Pwave'*Pwave*Vh;      %ArBC c related part 
    end
    
    
    
    
    
    
    
    %% form B vector basis     
    
    %Br=Brk*Zk;
    Brk=Vh'*P'*P * -(UpShift1Eye-lowShift1Eye)./(2*deltaZ)*VdK;   
    %This structure allows POD and DEIM have different number of basis.
    
    for i=1:size(Vh,2)  %number of pod basis
       Brhc(:,:,i)=Vh'*P'*P*spdiags(Vh(:,i),0,nZ,nZ)*VdC;       
%        Br         = Br +Brhc(:,:,i)*previousZh(i)*Zc./deltaT;
    end
    
    %Original method would slow down pod model thus not is not used
%     Br =   Vh'*P'*P * -(UpShift1Eye-lowShift1Eye)./(2*deltaZ)*Vk*Zk...
%           +Vh'*P'*P *   ((Vh*previousZh) ./deltaT .*(Vc *Zc));
      
    %%Take away BC point in ROM  
    %Br =   Br- Ar2*Zh;
    
    %TODO the following structure could be improve to save memory storage.
    %% Save rom model
    romMesh.Ark=Ark;
    romMesh.Arc=Arc;

    romMesh.ArBCk=ArBCk;
    romMesh.ArBCc=ArBCc;

    romMesh.Brk=Brk;
    romMesh.Brhc=Brhc;
    
    romMesh.nPod=nPod;
    romMesh.nDeimK=nDeimK;
    romMesh.nDeimC=nDeimC;
    
    % another way to accelerate
    romMesh.mArk=reshape(Ark,[],nDeimK);
    romMesh.mArc=reshape(Arc,[],nDeimC);

    romMesh.mArBCk=reshape(ArBCk,[],nDeimK);
    romMesh.mArBCc=reshape(ArBCc,[],nDeimC);

    romMesh.mBrk=reshape(Brk,[],nDeimK);
    romMesh.mBrhc=reshape(Brhc,[],nPod);
    
    romMesh.Vh=Vh;
    romMesh.Zh=Vh'*mesh.H;
    romMesh.Ks=mesh.Ks;
    
    romMesh.Pk=Pk;      %this oculd be easily compressed
    romMesh.Pc=Pc;  
    
%     romMesh.Brhc2=Brhc;
    
    
%         Ar=Ar+Ark(:,:,i).*Zk(i) +Arc(:,:,i).*Zc(i);
% 
%         Br=Brk*Zk;
%         Br         = Br +Brhc(:,:,i)*previousZh(i)*Zc;
% 
%         ArBC=Ar+ArBCk(:,:,i).*Zk(i) +ArBCc(:,:,i).*Zc(i);
% 
%         Brk=Vh'*P'*P * -(UpShift1Eye-lowShift1Eye)./(2*deltaZ)*Vk;
% 
%         %last modify for BC
%         Br =   Br- Ar2*Zh;
end


