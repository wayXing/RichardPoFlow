function [H,iteration] = Richard3dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K)
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
% History:  01/06/2017  file created
% TODO:     -Accept non-uniform t input
%           -only accept homogeneous Neumann BC or Dirichlet BC.


%% Auxiliary variable
nodeInFieldIndex=find(mesh.nodeIndex); %specify free node index  
%% Main
for iT=1:nTime
    previousH=mesh.H;
      % Picard iteration
    for k=1:nMaxIteration 
        H0=mesh.H;  %preserved for iteration compare
        
        %update mesh value 
%         mesh.C=theataDifFunc(mesh.H);
% %         mesh.K=kFunc(mesh.H);
%         mesh.K=kFieldFunc(mesh.H,mesh.Ks);
        
        mesh.C=theataDif(mesh.H);
        mesh.K=K(mesh.H,mesh.Ks);
        
        [A,B]=picard3dAxbForm2(mesh,previousH,deltaT);
        %solve linear equation
        h=A\(-B);
        
        %update mesh value
        mesh.H(nodeInFieldIndex)=h; 
         
%         hSnapShot=[hSnapShot,mesh.H];   %maybe used as snapshots

        %stopping criteria
        sseIte=sum((mesh.H(:)-H0(:)).^2);
        if sqrt(sseIte)<maxIteError 
            iteration(iT,1)=k;     %record n iteration till converge
            break 
        else 
%             nIterationFom=nIterationFom+1;
        end
    end 
    H(:,:,:,iT)=mesh.H;
    
end




end