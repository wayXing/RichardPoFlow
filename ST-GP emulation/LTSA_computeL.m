function [ model ] = LTSA_computeL( model )
% Local tangent space alignment preimage
% 
% Description:
%     This function is computes the local affine transformation matrices L
%     and their respective inverse for each training data point.
%     
%   Input:     
%       model [struct] Decribes the parameter of LTSA
%             .DR_method [string]       the method of reducing dimension
%             .options                        keep the input data
%             .cputime                        process time
%             .X                              original data
%             .k                              number of neighbor point
%             .S                              slection matrix
%             .theta_inv                      inverse matrix of theta
%             .ni                             neighbor point list
%             .Q                              d left singular vector    
%   Output:
%       model [struct] Decribes the parameter of LTSA
%             .DR_method [string]       the method of reducing dimension
%             .options                        keep the input data
%             .cputime                        process time
%             .X                              original data
%             .k                              number of neighbor point
%             .S                              slection matrix
%             .theta_inv                      inverse matrix of theta
%             .ni                             neighbor point list
%             .Q                              d left singular vector    
%             .L                              local affine transformation matrices
%             .L_inv                           inv(L)
% 
% See also
%     LTSA
% 
% About:
%     Modification
%     Zheng Xing,11-8-2016,First Edition
% 
% Reference:
% PRINCIPAL MANIFOLDS AND NONLINEAR DIMENSION REDUCTION VIA LOCAL TANGENT 
% SPACE ALIGNMENT 2004,written by Zhenyue Zhang & Hongyuan Zha

%% Initialization and Parameters
    k=model.neighbor;
    S=model.S;
    
    T=model.T;
    T=real(T);              % should probably have asserted this rather than convert it when it has already failed...
    assert( isreal(T) );    % added assertion to training so this should be redundant
   

%% main

% compute the inverse of L,refer to L- 
    % pre-alloc
    L     = zeros(size(model.X,1), size(model.T,2), size(model.T,2));
    L_inv =  zeros(size(model.X,1), size(model.T,2), size(model.T,2));
    Ti = zeros(size(model.X,1),size(model.T,2),k);

    h = waitbar(0,'Calculating matrix inverses...');
    for i=1:size(model.X,1)
        waitbar(i/size(model.X,1), h);
        Ti(i,:,:) = T'*S(:,(k*(i-1)+1):(k*i));
        theta_inv_temp = squeeze( model.theta_inv(i,:,:) );        
        L(i,:,:)= squeeze(Ti(i,:,:))*(eye(k)-ones(k,k)./k)*theta_inv_temp;
        L_inv(i,:,:)=inv( squeeze(L(i,:,:)) );
    end

    model.Ti = Ti;
    model.L = L;
    model.L_inv = L_inv;


