function [A,B]=picard3dAxbForm2(mesh,H_PreviousTime,deltaT)
%Initilize Picard iteration using DEIM ROM
%on the 1d h-based Richards equation
%
% Input parameters:
%   mesh             -mseh structure 
% Output parameters:
%   Ar,Br            -A*Z=B;
%
% Examples: see Demo
%
% See also: 
% Author:   Wei Xing
% History:  31/05/2017  file created
% log:      
% v2        -created from picard3dAxbForm to speed up using better storage
%           structure (avoid visiting full, though sparse, matrix at every 
%           iteration. NO proper vectorization yet. Loop in use!

%%  Auxiliary variable   
nNode=mesh.nNode;
deltaZ=mesh.deltaZ;
nZ=mesh.nZ;

deltaX=mesh.deltaX;
nX=mesh.nX;

deltaY=mesh.deltaY;
nY=mesh.nY;

nodeIndex=mesh.nodeIndex;


H=mesh.H;
C=mesh.C;
K=mesh.K;

% C=theataDifFunc(H);
% K=kFieldFunc(H,Ks);
A=speye(length(nNode));
B=zeros(length(nNode),1);


%     for iZ=1:nZ
%         for iX=1:nX
%             for iY=1:nY
                
   
    
    for iY=1:nY
        for iX=1:nX
            for iZ=1:nZ
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

%                 cCenter=C(iZ,1); 
                cCenter=C(iZ,iX,iY);

                wUp   = -kHalfUp  ./deltaZ^2;
                wDown = -kHalfDown./deltaZ^2;
                wLeft = -kHalfLeft./deltaX^2;
                wRight= -kHalfRight./deltaX^2;                    
                wFront= -kHalfFront./deltaY^2;
                wBack = -kHalfBack ./deltaY^2;


                wCenter=cCenter/deltaT-wUp-wDown-wLeft-wRight-wFront-wBack;

                b=(kHalfDown-kHalfUp)/deltaZ-H_PreviousTime(iZ,iX,iY)*cCenter/deltaT;

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
                
                iMethod=2;
                switch iMethod
                    case 1  %slow method as A is large.
                        if indexUp>0 A(indexCenter,indexUp)=wUp; end
                        if indexDown>0 A(indexCenter,indexDown)=wDown; end
                        if indexLeft>0 A(indexCenter,indexLeft)=wLeft; end
                        if indexRight>0 A(indexCenter,indexRight)=wRight; end
                        if indexFront>0 A(indexCenter,indexFront)=wFront; end
                        if indexBack>0 A(indexCenter,indexBack)=wBack; end                    

                        A(indexCenter,indexCenter)=wCenter;
                        B(indexCenter,1)=b;
                    case 2
                        upBand(indexCenter,1)    =wUp    *(indexUp>0);
                        downBand(indexCenter,1)  =wDown  *(indexDown>0);
                        leftBand(indexCenter,1)  =wLeft  *(indexLeft>0);
                        rightBand(indexCenter,1) =wRight *(indexRight>0);
                        frontBand(indexCenter,1) =wFront *(indexFront>0);
                        backBand(indexCenter,1)  =wBack  *(indexBack>0);
                        centerBand(indexCenter,1)=wCenter;
                        
                        B(indexCenter,1)=b;
                end   
            end %iZ loop
        end %iX loop
    end %iY loop
        
    
    %use when using iMethod 2 NO need for iMethod 1
    %TODO: need structure(mesh) refinement to improve. Only work on special
    %BC here. (top&bottom DBC and sides NBC)
    nFreeZ=nZ-2;
    nFreeX=nX;
    nFreeY=nY;
    
    A=spdiags([downBand,upBand,rightBand,leftBand,backBand,frontBand,centerBand],[-1,1,-nFreeZ,nFreeZ,-nFreeZ*nFreeX,nFreeZ*nFreeX,0],nNode,nNode);
    
end

