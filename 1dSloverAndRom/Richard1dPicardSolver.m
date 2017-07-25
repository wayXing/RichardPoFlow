function [H,iteration] = Richard1dPicardSolver(mesh,nTime,deltaT,nMaxIteration,maxIteError,theataDif,K)
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
nodeIndex=find(~mesh.dbcFlag);   %specify free node index   

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
        
        [A,B]=picardAxbForm(mesh,previousH,deltaT);
        %solve linear equation
        h=A\(B);
        
        %update mesh value
        mesh.H(nodeIndex)=h;
         
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
    H(:,iT)=mesh.H;
    
end




end