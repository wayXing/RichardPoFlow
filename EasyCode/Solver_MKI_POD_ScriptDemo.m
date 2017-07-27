function [] = Solver_MKI_POD_ScriptDemo()
% 1D Richards equation with constant boundary condition 
%
% This is a self-contained Script DEMO function that requires no external
% functions and data.
% This Code is meant for tutorial purpose and thus may suffer from
% computational inefficiency and other issues.
%
% First edition: Weix 23/03/2017 
%%
clear 

nNode=10;
nTime=100;

% Spatial setup
lengthZ=40;
deltaZ=1;
nNode=lengthZ/deltaZ+1;

% Temporal setup
lengthTime=360;
deltaTime=1;
nTime=lengthTime/deltaTime+1;

% Iteration setup
nMaxIteration=1000;
miniIteError=0.01;



H=zeros(nNode,nTime);
H_temp=zeros(nNode,nMaxIteration);

%%% 

% H=zeros(nTime,nNode,nMaxIteration);
% time=[0:deltaTime:100]

% Initial condition and Boundary condition
h_init=zeros(nNode,1);
h_init(:)=-61.5;
h_init(1)=-20.7;
h_init(end)=-61.5;


%% MAIN
tic

H(:,1)=h_init;
h=h_init;
H_all=[];

for t=1:nTime
    
    h_PreviousTime=h;
    h0=h_PreviousTime; %???
    
    for j=1:nMaxIteration     

        %Assemble matrix for A_l*x=b_l
        matrixA=zeros(nNode);
        columnB=zeros(nNode,1);
        
        for i = 2:nNode-1
                
            %c by SCS approximation
%             h0(i)=h0(i)+0.0000001;    %To avoid NaN      
            c(i)=(theata(h0(i))-theata(h_PreviousTime(i)))/(h0(i)+0.00000001-h_PreviousTime(i)); % add 0.0000001 To avoid NaN
            
            %c by analytical derivatives
%             c(i)=theataDif(h(i));
            
            
            a=(k(h0(i))+k(h0(i-1)))/(-2*deltaZ^2);
            b=c(i)/deltaTime+(k(h0(i+1))+2*k(h0(i))+k(h0(i-1)))/(2*deltaZ^2);
            d=-(k(h0(i+1))+k(h0(i)))/(2*deltaZ^2);
            e=(-k(h0(i+1))+k(h0(i-1)))/(2*deltaZ)+c(i)*h_PreviousTime(i)/deltaTime;

            matrixA(i,i-1)=a;
            matrixA(i,i)=b;
            matrixA(i,i+1)=d;

            columnB(i)=e;

        end
        
        %Deal with B.C.(truncate the system matrix and take care of the
        %consequence)
        
        columnB=columnB-h_init(1).*matrixA(:,1);
        columnB=columnB-h_init(end).*matrixA(:,end);
        
        
        matrixA=matrixA(2:nNode-1,2:nNode-1);
        columnB=columnB(2:nNode-1,1);
        
%         columnB(1)=columnB(1)-matrixA(1,2)*h_init(1);
%         columnB(end)=columnB(end)-matrixA(end-1,end)*h_init(end);
        
       
        %Utilize sparsity
%         sparseA=sparse(matrixA);
%         h=sparseA\columnB;
        
        
        %Solve linear system
        h=matrixA\columnB;
        
        %Complete h with B.C.
        h=[h_init(1);h;h_init(end)];
        
        
        smeIte=sum((h-h0).^2);
        if sqrt(smeIte)<miniIteError 
            break 
        else
            h0=h;
        end
        
        H_all=[H_all,h];

    end
    
    H(:,t+1)=h;
    
end
computeTime=toc

%% POD 
% nPOD=10;
podEnergy=0.999;

%better use all 
[U,S,V]=svd(H_all);
% [U,S,V]=svd(H);
energy=diag(S);
cumulatedEnergy= cumsum(energy)./sum(energy);
[~,nPOD]=min(abs(cumulatedEnergy-podEnergy))

podBasis=U(2:end-1,1:nPOD);


% %better use all 
% [U,S,V]=svd(H_all);
% % [U,S,V]=svd(H);
% 
% podBasis=U(2:end-1,1:nPOD);
    

H2(:,1)=h_init;
h=h_init;

tic 
for t=1:nTime
    
    h_PreviousTime=h;
    h0=h_PreviousTime; %???
    
    for j=1:nMaxIteration     

        %Assemble matrix for A_l*x=b_l
        matrixA=zeros(nNode);
        columnB=zeros(nNode,1);
        
        for i = 2:nNode-1
                
            %c by SCS approximation
%             h0(i)=h0(i)+0.0000001;    %To avoid NaN      
            c(i)=(theata(h0(i))-theata(h_PreviousTime(i)))/(h0(i)+0.00000001-h_PreviousTime(i)); % add 0.0000001 To avoid NaN
            
            %c by analytical derivatives
%             c(i)=theataDif(h(i));
            
            
            a=(k(h0(i))+k(h0(i-1)))/(-2*deltaZ^2);
            b=c(i)/deltaTime+(k(h0(i+1))+2*k(h0(i))+k(h0(i-1)))/(2*deltaZ^2);
            d=-(k(h0(i+1))+k(h0(i)))/(2*deltaZ^2);
            e=(-k(h0(i+1))+k(h0(i-1)))/(2*deltaZ)+c(i)*h_PreviousTime(i)/deltaTime;

            matrixA(i,i-1)=a;
            matrixA(i,i)=b;
            matrixA(i,i+1)=d;

            columnB(i)=e;
            
            

        end

        
 
        
        %Deal with B.C.(truncate the system matrix and take care of the
        %consequence)
        
        columnB=columnB-h_init(1).*matrixA(:,1);
        columnB=columnB-h_init(end).*matrixA(:,end);
        
        
        matrixA=matrixA(2:nNode-1,2:nNode-1);
        columnB=columnB(2:nNode-1,1);
        
%         columnB(1)=columnB(1)-matrixA(1,2)*h_init(1);
%         columnB(end)=columnB(end)-matrixA(end-1,end)*h_init(end);
        
       
        %Utilize sparsity
%         sparseA=sparse(matrixA);
%         h=sparseA\columnB;
        
        
        %MOR 
        matrixAr=podBasis'*matrixA*podBasis;
        columnBr=podBasis'*columnB;




        %Solve linear system
        hr=matrixAr\columnBr;
        
        h=podBasis*hr;
        %Complete h with B.C.
        h=[h_init(1);h;h_init(end)];
        
        
        smeIte=sum((h-h0).^2);
        if sqrt(smeIte)<miniIteError 
            break 
        else
            h0=h;
        end
        

    end
    
    H2(:,t+1)=h;
    
end
podComputeTime=toc

figure(1)
for t=1:nTime
    plot(H(:,t))
    hold on
    plot(H2(:,t))
    hold off
    
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;
    
end







    
end
    
    

function result=theata(h)
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

result=alpha*(theata_s-theata_r)/(alpha+abs(h)^beta)+theata_r;
end

function result=theataDif(h)
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

result=-alpha*(theata_s-theata_r)*-1*(alpha+abs(h)^beta)^(-2)*abs(h)^(beta-1);

end

function result=k(h)
rho=1.175e6;
r=4.74;
k_s=0.00944;

result=k_s*rho/(rho+abs(h)^r);
end



