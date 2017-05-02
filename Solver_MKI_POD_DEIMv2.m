function [H,H_POD] = Solver_MKI_POD_DEIMv2()
% 1D Richards equation with constant boundary condition 
%
% First edition: Weix 23/03/2017 
% v2: Weix 23/03/2017 
% improve speed. using permibility field. use podEnergy than nPod.

% nNode=10;
% nTime=100;

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

a2_record=[];
b_record=[];
d_record=[];
e_record=[];

for t=1:nTime
    
    h_PreviousTime=h;
    h0=h_PreviousTime; %???
    
    for j=1:nMaxIteration     

        %Assemble matrix for A_l*x=b_l
        matrixA=zeros(nNode);
        columnB=zeros(nNode,1);
        
        a=zeros(nNode-1,1);
        a2=zeros(nNode-1,1);
        
        b=zeros(nNode,1);
        d=zeros(nNode-1,1);
        e=zeros(nNode-2,1);
        
        %Assemble the A matrix and b column
        for i = 2:nNode-1
                
            %c by SCS approximation
%             h0(i)=h0(i)+0.0000001;    %To avoid NaN      
            c(i)=(theata(h0(i))-theata(h_PreviousTime(i)))/(h0(i)+0.00000001-h_PreviousTime(i)); % add 0.0000001 To avoid NaN
            
            %c by analytical derivatives
%             c(i)=theataDif(h(i));
            
            
            a(i)=(k(h0(i))+k(h0(i-1)))/(-2*deltaZ^2);
            a2(i-1)=(k(h0(i))+k(h0(i-1)))/(-2*deltaZ^2);
            
            b(i)=c(i)/deltaTime+(k(h0(i+1))+2*k(h0(i))+k(h0(i-1)))/(2*deltaZ^2);
            d(i)=-(k(h0(i+1))+k(h0(i)))/(2*deltaZ^2);
            e(i)=(-k(h0(i+1))+k(h0(i-1)))/(2*deltaZ)+c(i)*h_PreviousTime(i)/deltaTime;

            % One way to assemble
            matrixA(i,i-1)=a(i);
            matrixA(i,i)=b(i);
            matrixA(i,i+1)=d(i);

            columnB(i)=e(i);

        end
        
        a2_record=[a2_record,a2];
        b_record=[b_record,b];
        d_record=[d_record,d];
        e_record=[e_record,e];
        
        
        %Another way to assemble
        matrixA2=diag(a2,-1)+diag(b,0)+diag(d,1);
        
%         sum(matrixA2(:)-matrixA(:))
        matrixA=matrixA2;
        
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

%% DEIM pre-compute 
nDeim=30;

[Ua,~,~]=svd(a2_record);
[~,~,Pa] = DEIM(Ua);
Pra=Pa(:,1:nDeim);
Ura=Ua(:,1:nDeim);
Dra=Ura*inv(Pra'*Ura);

[Ub,~,~]=svd(b_record);
[~,~,Pb] = DEIM(Ub);
Prb=Pb(:,1:nDeim);
Urb=Ub(:,1:nDeim);
Drb=Urb*inv(Prb'*Urb);

[Ud,~,~]=svd(d_record);
[~,~,Pd] = DEIM(Ud);
Prd=Pd(:,1:nDeim);
Urd=Ud(:,1:nDeim);
Drd=Urd*inv(Prd'*Urd);






%% POD 
nPOD=10;

a2r_rec=[];
a2_rec=[];

%better use all 
[U,S,V]=svd(H_all);
% [U,S,V]=svd(H);

podBasis=U(2:end-1,1:nPOD);
    

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
        
        a=zeros(nNode-1,1);
        a2=zeros(nNode-1,1);
        
        b=zeros(nNode,1);
        d=zeros(nNode-1,1);
        e=zeros(nNode-2,1);
        
        %Assemble the A matrix and b column
        for i = 2:nNode-1
                
            %c by SCS approximation
%             h0(i)=h0(i)+0.0000001;    %To avoid NaN      
            c(i)=(theata(h0(i))-theata(h_PreviousTime(i)))/(h0(i)+0.00000001-h_PreviousTime(i)); % add 0.0000001 To avoid NaN
            
            %c by analytical derivatives
%             c(i)=theataDif(h(i));
            
            
            a(i)=(k(h0(i))+k(h0(i-1)))/(-2*deltaZ^2);
            a2(i-1)=(k(h0(i))+k(h0(i-1)))/(-2*deltaZ^2);
            
            b(i)=c(i)/deltaTime+(k(h0(i+1))+2*k(h0(i))+k(h0(i-1)))/(2*deltaZ^2);
            d(i)=-(k(h0(i+1))+k(h0(i)))/(2*deltaZ^2);
            e(i)=(-k(h0(i+1))+k(h0(i-1)))/(2*deltaZ)+c(i)*h_PreviousTime(i)/deltaTime;

            % One way to assemble
            matrixA(i,i-1)=a(i);
            matrixA(i,i)=b(i);
            matrixA(i,i+1)=d(i);

            columnB(i)=e(i);

        end
        
        a2r= Dra*(Pra'*a2);
        br= Drb*(Prb'*b);
        dr= Drd*(Prd'*d);
        
        a2r_rec=[a2r_rec,a2r];
        a2_rec=[a2_rec,a2];
        
        
        
%         er= Dre*(Pre'*e);
        
        
%         sum(a2r(:)-a2(:))

%         figure(1)
%         plot(a2r)
%         hold on
%         plot(a2)
%         
%         plot(br)
%         plot(b)
%         
%         plot(dr)
%         plot(d)
%         
%         hold off
        
        %Another way to assemble
        matrixA2=diag(a2r,-1)+diag(br,0)+diag(dr,1);
%         sum(matrixA2(:)-matrixA(:))
        matrixA=matrixA2;
        

        

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



% figure(1)
% for i=1:1000
%     plot(a2r_rec(:,i))
%     hold on
%     plot(a2_rec(:,i))
%     hold off
%     drawnow
%     frame(i)=getframe;
%     
% end






figure(2)
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