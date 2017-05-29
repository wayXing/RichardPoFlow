function [Ar,Br]=picardAxbRomForm(romMesh,deltaT,previousZh,Zh,Zk,Zc)
%Form linear system at picards iteration on a initilized reduced order model 
%using DEIM ROM on the 1d h-based Richards equation
%
% Input parameters:
%   romMesh             -initilized reduced order model
%   deltaT              -
%   previousZh          -value of Zh at last time step
%   Zh                  -current Zh
%   Zk                  -current Zk
%   Zc                  -current Zc
%
% Output parameters:
%   Ar,Br            -Ar*Zh=Br;
%
% Examples: see Demo
%
% See also: 
% Author:   Wei Xing
% History:  22/05/2017  file created

nPod=romMesh.nPod;
nDeimK=romMesh.nDeimK;
nDeimC=romMesh.nDeimC;


Ar=zeros(nPod,nPod);
ArBC=zeros(nPod,nPod);

%Ar matrix 
% for i=1:nDeimK
% 
%     Ar= Ar+romMesh.Ark(:,:,i).*Zk(i) +romMesh.Arc(:,:,i).*Zc(i)./deltaT;
%     
%     ArBC= ArBC+romMesh.ArBCk(:,:,i).*Zk(i) +romMesh.ArBCc(:,:,i).*Zc(i)./deltaT;        
% 
% end

iMethod=2;
switch iMethod 
    case 1  %Normal process
        %Ar matrix k related term
        for i=1:nDeimK
            Ar  = Ar   +romMesh.Ark(:,:,i).*Zk(i);
        %     Ar  = Ar   + bsxfun(@times, romMesh.Ark(:,:,i),Zk(i));   
            ArBC= ArBC +romMesh.ArBCk(:,:,i).*Zk(i);       
        end

        %Ar matrix c related term
        for i=1:nDeimC
            Ar  = Ar   +romMesh.Arc(:,:,i).*Zc(i)./deltaT;
            ArBC= ArBC +romMesh.ArBCc(:,:,i).*Zc(i)./deltaT;        
        end
        
        %Improve calculation by using gpu parellel
%         tempZc=reshape(Zc./deltaT,1,1,length(Zc));
%         Ar=sum(pagefun(@mtimes,romMesh.Brhc,tempZc),3);

        %Improve calculation by not using for loop
        mArk=reshape(romMesh.Ark,nDeimK*nDeimK,nDeimK);
        mArc=reshape(romMesh.Arc,nDeimK*nDeimK,nDeimK);
        Ar=mArk*Zk+mArc*Zc;
        Ar=reshape(Ar,nDeimK,nDeimK);

        %Br matrix
        Br=romMesh.Brk*Zk;
        for i=1:nPod  %number of pod basis
            Br = Br +romMesh.Brhc(:,:,i)*(previousZh(i)./deltaT*Zc);
        end

%         tempZh=reshape(previousZh,1,1,length(previousZh));
%         temp=previousZh./deltaT*Zc;
%         temp=reshape(temp,[],1,length(previousZh));
%         
%         Br=romMesh.Brk*Zk;
% %         Br = Br +MULTIPROD(romMesh.Brhc, (previousZh./deltaT*Zc'));
%         Br = Br +pagefun(@times,romMesh.Brhc,(tempZh./deltaT*Zc));
        
        %%Take away BC point in ROM  
        Br =   Br- ArBC*Zh;
        
        
    case 2 %vectorize to speed up
        Ar=romMesh.mArk*Zk+romMesh.mArc*Zc./deltaT;
        ArBC=romMesh.mArBCk*Zk+romMesh.mArBCc*Zc./deltaT;
        
        Ar=reshape(Ar,nPod,nPod);
        ArBC=reshape(ArBC,nPod,nPod);
        
        %Br matrix
%         Br=romMesh.Brk*Zk;
%         for i=1:nPod  %number of pod basis
%             Br = Br +romMesh.Brhc(:,:,i)*(previousZh(i)./deltaT*Zc);
%         end
        
        Br2= romMesh.mBrhc*previousZh;
%         Br2= reshape(Br2,nPod,nPod);
        Br2= reshape(Br2,nPod,nDeimC);
        Br= Br2*Zc+romMesh.Brk*Zk;

        %%Take away BC point in ROM  
        Br =   Br- ArBC*Zh;
        
    case 3  %parellel process 
        %Ar matrix k related term
        parfor i=1:nDeimK
            Ar  = Ar   +romMesh.Ark(:,:,i).*Zk(i);
        %     Ar  = Ar   + bsxfun(@times, romMesh.Ark(:,:,i),Zk(i));   
            ArBC= ArBC +romMesh.ArBCk(:,:,i).*Zk(i);       
        end

        %Ar matrix c related term
        parfor i=1:nDeimC
            Ar  = Ar   +romMesh.Arc(:,:,i).*Zc(i)./deltaT;
            ArBC= ArBC +romMesh.ArBCc(:,:,i).*Zc(i)./deltaT;        
        end


        %Br matrix
        Br=romMesh.Brk*Zk;
        parfor i=1:nPod  %number of pod basis
            Br = Br +romMesh.Brhc(:,:,i)*(previousZh(i)./deltaT*Zc);
        end

        %%Take away BC point in ROM  
        Br =   Br- ArBC*Zh;
        
end




end