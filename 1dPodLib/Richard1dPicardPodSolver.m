function [H,iteration] = Richard1dPicardPodSolver(romMesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K,Pk,Pc)
% solve Richards equation in 1d with picard iteration.
% the Boundary condition should be included in the Mesh.
%
% Input parameters:
%
% Output parameters:
%
% See also: 
%
% Author:   Wei Xing
% History:  26/05/2017  file created
% TODO:     Accept non-uniform t input

%% Auxiliary variable
Vh=romMesh.Vh;

%% Initialize
iteration=ones(nTime,1).*nMaxIteration;


%% Main

VhTVh=Vh'*Vh;
Zh=romMesh.Zh;

% Pk=romMesh.Pk;
% Pc=romMesh.Pc;

PkVh=Pk'*Vh;
PcVh=Pc'*Vh;

tic
for iT=1:nTime   
    previousZh=Zh;
    % Picard iteration 
    for k=1:nMaxIteration 
        Zh0=Zh;  %preserved for iteration compare
         
        %update mesh value        
     
        %TODO simple this to speed up
%         Zk=K(Pk'*Vh*Zh,Pk'*romMesh.Ks);
%         Zc=theataDif(Pc'*Vh*Zh);
                       
        Zk=K(PkVh*Zh,Pk'*romMesh.Ks);
        Zc=theataDif(PcVh*Zh);        
        
        
        [Ar,Br]=picardAxbRomForm(romMesh,deltaT,previousZh,Zh,Zk,Zc);
             
        %solve linear equation
        Zh=Ar\(Br);
              
        %stopping criteria
%         mesh.H=V*Zh;
%         sseIte=sum((mesh.H(:)-H0(:)).^2);
        sseIte=(Zh-Zh0)'*VhTVh*(Zh-Zh0); 
        if sqrt(sseIte)<maxIteError 
            iteration(iT,1)=k;     %record n iteration till converge
            break 
        else 
%             nIterationDeimPod=nIterationDeimPod+1;
        end
%         disp('.')
    end 
    ZhRecord(:,iT)=Zh;
    
end
timeCost=toc;
%after process
H=Vh*ZhRecord;





end