function [ Y_star,model_2 ] = LTSA_preimage3( Tt,model )
% Local tangent space alignment preimage
% 
% Description:
%     This function is implementation of LTSA preimage.The input data is reduced
%     dimension by LTSA.The result is the reconsruction of X.
%     
%   Input:
%    	Tt [Num_data x dim] Training data.
%     
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
%       X_star [Num_data x org_dim]           The reconstruct of X
%     
%       model 
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
    [Num_data_test, new_dim]=size(Tt);
    k=model.neighbor;
    S=model.S;
    ni=model.ni;
    Q=model.Q;
    Ti = model.Ti;
    L_inv = model.L_inv;
    L = model.L;
    
    X_fix_ner_avg = model.X_fix_ner_avg;
    Ti_fix_ner_avg = model.Ti_fix_ner_avg;
    
    T=model.T;
    T=real(T);              % should probably have asserted this rather than convert it when it has already failed...
    assert( isreal(T) );    % added assertion to training so this should be redundant
    [Num_data_train, org_dim]=size(X);
    
    [ni2, ~]=knnsearch(T,Tt,'k',1,'distance','euclidean');
   
    model_2.ni2=ni2;

%% main

    % predict
    % pre-allo
    Y_star = zeros(Num_data_test, org_dim);
    for i=1:Num_data_test
        Y_star(i,:)=X_fix_ner_avg(ni2(i),:)+ ...
                    ( squeeze(Q(ni2(i),:,:))* squeeze(L_inv(ni2(i),:,:))*...
                      (Tt(i,:) - Ti_fix_ner_avg(ni2(i),:))'  )';
    end

