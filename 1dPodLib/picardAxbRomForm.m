function [Ar,Br]=picardAxbRomForm(romMesh,deltaT,previousZh,Zh,Zk,Zc)
%Initilize Picard iteration using DEIM ROM
%on the 1d h-based Richards equation
%
% Input parameters:
%   mesh             -mseh structure 
%   previousH        -value of mesh.H at last time step
%   Vh                  
%
% Output parameters:
%   A,B              -A*mesh.H_new=B;
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

%Ar matrix k related term
for i=1:nDeimK
    Ar  = Ar   +romMesh.Ark(:,:,i).*Zk(i);
    ArBC= ArBC +romMesh.ArBCk(:,:,i).*Zk(i);        
end

%Ar matrix c related term
for i=1:nDeimC
    Ar  = Ar   +romMesh.Arc(:,:,i).*Zc(i)./deltaT;
    ArBC= ArBC +romMesh.ArBCc(:,:,i).*Zc(i)./deltaT;        
end





%Br matrix
Br=romMesh.Brk*Zk;
for i=1:nPod  %number of pod basis
    Br = Br +romMesh.Brhc(:,:,i)*previousZh(i)*Zc./deltaT;
end

%%Take away BC point in ROM  
Br =   Br- ArBC*Zh;




end