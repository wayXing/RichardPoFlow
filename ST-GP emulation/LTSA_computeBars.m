function [ model ] = LTSA_computeBars( model )
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
    X=model.X;
    k=model.neighbor;
    ni=model.ni;
    Ti = model.Ti;
    T=model.T;
    T=real(T);              % should probably have asserted this rather than convert it when it has already failed...
    assert( isreal(T) );    % added assertion to training so this should be redundant
    
    new_dim = model.new_dim;
    [Num_data_train, org_dim]=size(X);



%% main

% compute Ti,refer to the Ti in paper
% pre-alloc

X_fix_ner_avg=zeros(Num_data_train,org_dim);
Ti_fix_ner_avg=zeros(Num_data_train,new_dim);
for i=1:Num_data_train
    for j=1:k
        X_fix_ner_avg(i,:)=X(ni(i,j),:)+X_fix_ner_avg(i,:);
        Ti_fix_ner_avg(i,:)=Ti(i,:,j)+Ti_fix_ner_avg(i,:);
    end
end
% compute X_bar and t_bar (center of tangent space)
X_fix_ner_avg=X_fix_ner_avg/k;              %the average of x's neighborhood (inc itself)
Ti_fix_ner_avg=Ti_fix_ner_avg/k;            %the average of t's neighborhood (inc itself)


model.X_fix_ner_avg  = X_fix_ner_avg;
model.Ti_fix_ner_avg = Ti_fix_ner_avg;


