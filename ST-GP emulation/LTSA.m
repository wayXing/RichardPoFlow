function [ T,model ] = LTSA( X,options )
%LTSA Local Tangent Space Alignment.
%     
% Description:
%     This Function is implementation of Local Tangent Space Alignment(LTSA).The input 
%     data X are non-linearly mapped to a new local tangent space.The result is a model
%     describing non-linear data projection.
% 
% To do:
%     1. Implement local summing for computation of B
%     2. Implement incremental manifold learning for ALM
%
% Input:
%     X [Num_data x dim] Training data.
%     
%     option [struct] Decribes 
%         .k [num]                            the number of neighbor point
%         .new_dim [num]                      the reduced dimension
%     Output:
%         T [Num_data x new_dim]              matrix of projection of X
%         model
%             .DR_method [string]       the method of reducing dimension
%             .options                        keep the input data
%             .cputime                        process time
%             .X                              original data
%             .k                              number of neighbor point
%             .S                              slection matrix
%             .theta_inv                      inverse matrix of theta
%             .ni                             neighbor point list
%             .Q                              d left singular vector
%             
% see also
%     find_nn
%   
% About:
%     Modification
%     Zheng Xing, 8-11-2016,First Edition
% 
% Reference:
% PRINCIPAL MANIFOLDS AND NONLINEAR DIMENSION REDUCTION VIA LOCAL TANGENT 
% SPACE ALIGNMENT 2004,written by Zhenyue Zhang & Hongyuan Zha



%% Initialization and Parameters
% timer
start_time = cputime;

% process input arguments
%-----------------------------------
% if nargin < 2, options = [];end
% if ~isfield(options,'k'), options.k = 12; end
% if ~isfield(options,'new_dim'), options.new_dim = 2; end

k=options.neighbor;
[Num_data, org_dim]=size(X);
new_dim=options.new_dim;


%% Main

%find the k nearest point,see find_nn
% [ni] is the number list of the neighbor point
[ni, ~] = knnsearch(X,X,'k',k,'distance','euclidean');

% neighbor point inclouding the reference point itself 
% ni(:,2:k)=ni;
% ni(:,1)=1:Num_data;

% construct neighbor point matrix,refer to Xi[Xi1,Xi2~Xik]
    % pre-alloc
    Xk = zeros(Num_data, k, org_dim);
    for i=1:Num_data
        for j=1:k
            Xk(i,j,:)=X(ni(i,j),:);
        end
    end
    
    fprintf(1,'cpu time after constructing neughbour point matrix: %9.1f \n', cputime- start_time);

% centering the neighbor point matrix
% Xkc=X(I-ee'/k)
    % pre-allo
    Xkc = zeros(Num_data, org_dim, k);
    for i=1:Num_data
        Xkc(i,:,:)=squeeze(Xk(i,:,:))'*(eye(k)-ones(k,k)./k);
    end 
    
    fprintf(1,'cpu time after centering neighbour point matrix: %12.1f \n', cputime- start_time);
    
% compute Q,the d left singular vector of Xkc. refer to Qi in paper.
    Q = zeros(Num_data, org_dim, new_dim);
    for i=1:Num_data
        [U_temp, ~, ~] = svds( squeeze(Xkc(i,:,:)), new_dim );         %see svd function,S[eigen value]
        Q(i,:,:)=U_temp(:,1:new_dim);
    end  
    fprintf(1,'cpu time after computing Q: %33.1f \n', cputime - start_time);
    
    %compute theta and theta_inverse,refer to theta and theta+ 
    % pre-alloc
    theta = zeros(Num_data, new_dim, k);
    theta_inv = zeros(Num_data, k, new_dim);
    for i=1:Num_data    
        theta(i,:,:) = squeeze(Q(i,:,:))' * squeeze(Xkc(i,:,:));
        theta_inv(i,:,:)=pinv( squeeze(theta(i,:,:)) );
    end
    fprintf(1,'cpu time after computing theta and theta+: %18.1f \n', cputime- start_time);
    
% compute Wi
    for i=1:Num_data
        Wi(i,:,:)=(eye(k)-ones(k,1)*ones(1,k)./k)*(eye(k)-squeeze(theta_inv(i,:,:))* squeeze(theta(i,:,:))  );
    end
    fprintf(1,'cpu time after computing Wi: %32.1f \n', cputime- start_time);

% compute selection matrix S
    S=zeros(Num_data,k*Num_data);
    for i=1:Num_data
        S_temp=zeros(Num_data,k);
        for j=1:k
        S_temp(ni(i,j),j)=1;
        end
        S(1:Num_data,(k*(i-1)+1):(k*i))=S_temp;
    end
    fprintf(1,'cpu time after computing selection matrix: %18.1f \n', cputime- start_time);
 
% generalize W,diag with Wi   
    W=zeros(Num_data*k,Num_data*k);
    for i=1:Num_data
        W((k*(i-1)+1):k*i,(k*(i-1)+1):k*i)=Wi(i,:,:);
    end
    fprintf(1,'cpu time after concatinating W: %29.1f \n', cputime- start_time);
 
% compute B 
    B=S*W*W'*S';                                                                    % this needs decomposing
    fprintf(1,'cpu time after computing B: %33.1f \n', cputime- start_time);

    [E_vec, E_val]=eig(B);
    [~,ordered]=sort(diag(E_val));
    E_vec=E_vec(:,ordered);
    fprintf(1,'cpu time after finding T: %35.1f \n', cputime- start_time);
    
    T = E_vec(:,2:new_dim+1);
% [Tk_org,ni2]=find_nn(T,k-1);
% ni2(:,2:k)=ni2;
% ni2(:,1)=1:Num_data;


%% Recording & Output

model.DR_method='LTSA';
model.options = options;

model.cputime = cputime - start_time;
model.X=X;
model.neighbor=k;
model.S=S;
model.theta_inv=theta_inv;
model.ni=ni;
model.Q=Q;
model.B=B;
model.E_vec=E_vec;
model.W=W;
model.Num_data=Num_data;
model.org_dim=org_dim;
model.new_dim= new_dim;
model.theta=theta;
% model.ni2=ni2;
model.T=T;

model = LTSA_computeL(model);
model = LTSA_computeBars(model);

end

