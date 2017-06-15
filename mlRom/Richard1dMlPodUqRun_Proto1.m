function [] = Richard1dMlPodUqRun_Proto1()
% Initilize UQ for Richars equation with random input using mchaine learning & POD
% run training testing dataset
%
% Richars equation 1D pod solver testing file.
% The function focus on fix Dirichlet BC.
% This function serves as a Demo for all Richard solver developed in this
% project.
% 
% Proto1: Created form Richard1dMlPodUQ_Proto1() in idPod
%
% Input parameters:
%
% Output parameters:
%
% See also: 
%
% Author:   Wei Xing
% History:  15/06/2017  file created
%
clear
close all

filename='/Users/weix/Desktop/Data/data1';
% [status,msg] = mkdir('/uesrs/desktop/Data/data1.mat')
load(filename);

trainIndex=1:20;
testIndex =191:200;

HTrain=0;

HTrain=H_uq1(:,:,trainIndex);
xTrain=sample(:,trainIndex)';

HTest=H_uq1(:,:,testIndex);
xTest=sample(:,testIndex)';

%% Clustering
nCluster=3;

H_uq1Vec=reshape(HTrain,[],length(trainIndex));
labelTrian = kmeans(H_uq1Vec',nCluster);

%% local basis basend on Clusters
nPod=30;

h=waitbar(0,'Initilizing local ROMs');
for i=1:nCluster 
    
    index=find(labelTrian==i);
    
    % define snapshot
    iHSnapShot=reshape(HTrain(:,:,index),nZ,length(index)*nTime);   %decide snapshot
    
    %POD basis
%     [V_uq(:,:,i),S,~]=svds(iHSnapShot,nPod);
    [V_uq,S,~]=svds(iHSnapShot,nPod);
    
    % DEIM nonlinear function 
    nDeimK=nPod;    %number of Deim basis for k term
    nDeimC=nPod;    %number of Deim basis for c term
    
    %k
%     for t=1:size(iHSnapShot,2)
%         kRecord(:,t)=K(iHSnapShot(:,t),Ks(:,i));
%     end
    for j=1:length(index)
        for t=1:nTime
            kRecord(:,t,j)=K(HTrain(:,t,index(j)),Ks(:,index(j)));
        end  
    end
    kRecord=reshape(kRecord,nZ,[]);   %decide snapshot
    
    
    % [Vk,~,~]=svd(kRecord,'econ');
    [Vk,~,~]=svds(kRecord,nDeimK);

    [~,~,Pk] = DEIM(Vk);
    Pk=Pk(:,1:nDeimK);
    Vk=Vk(:,1:nDeimK);
%     VdK_uq(:,:,i)=Vk*inv(Pk'*Vk);  %DEIM basis
    VdK_uq=Vk*inv(Pk'*Vk);  %DEIM basis
    
    %c
    %     disp('DEIM decomposition for c...')
    cRecord=theataDif(iHSnapShot);
    % [Vc,~,~]=svd(cRecord,'econ');
    [Vc,~,~]=svds(cRecord,nDeimC);

    [~,~,Pc] = DEIM(Vc);
    Pc=Pc(:,1:nDeimC);
    Vc=Vc(:,1:nDeimC);
%     VdC_uq(:,:,i)=Vc*inv(Pc'*Vc);  %DEIM basis
    VdC_uq=Vc*inv(Pc'*Vc);  %DEIM basis
    
    romMesh{i}=picardAxbRomInit(mesh,V_uq,VdK_uq,Pk,VdC_uq,Pc);
    V_uqRecord(:,:,i)=V_uq;
    waitbar(i/nCluster)
end
close(h)

%% 
treeModel = fitctree(xTrain,labelTrian);
imp = predictorImportance(treeModel);
bar(imp);
title('Predictor Importance Estimates');
ylabel('Estimates');
xlabel('Predictors');

labelTest = predict(treeModel,xTest);

%% Deim POD
h=waitbar(0,'Deim pod on Ks on progress');
for i=1:length(testIndex)
    % Initilize ROM
    
    
    romMesh{labelTest(i)}.Ks=Ks(:,testIndex(i));
    
%     mesh.H=H_uq1(:,2,i); % use fom to start

%     romMesh{label(i)}.Zh=V_uqRecord(:,:,label(i))'*H_uq1(:,1,i);  %use FOM to start
    romMesh{labelTest(i)}.Zh=V_uqRecord(:,:,labelTest(i))'*h_init;
    
    tic
    [H_pod,iteration2] = Richard1dPicardPodSolver(romMesh{labelTest(i)},nTime,deltaT,nMaxIteration,maxIteError,theataDif,K);
    
    tCost2(i,1)=toc;
    iTera2(i,2)=sum(iteration2);
    H_uq2(:,:,i)=H_pod;
    
    waitbar(i/nSample)
end
close(h)


sum(tCost1)
sum(tCost2)

%% UQ process
mu_H_uq1 =mean(HTest,3);
var_H_uq1=std(HTest,0,3);
mid_H_uq1=median(HTest,3);

mu_H_uq2 =mean(H_uq2,3);
var_H_uq2=std(H_uq2,0,3);
mid_H_uq2=median(H_uq2,3);


%% show basis
% v1=reshape(V_uq(:,1,:),nZ,nSample);
% figure(1)
% plot(v1)



%% plot
nZShow=100;
zShow=1:round(nZ/nZShow):nZ;
figure(2)
for t=1:1:nTime
    figure(3)
    plot(squeeze( HTest(zShow,t,:)),'-')
    hold on 
    plot(squeeze( H_uq2(zShow,t,:)),'--')
    hold off
    ylim([-80,20])
    
    title(sprintf('time=%i',t))
%     legend('All KL basis','Truncation KL basis')
    drawnow
%     frame(t)=getframe;    %comment to save cpu time and memory
end


figure(3)
nZShow=100;
zShow=1:round(nZ/nZShow):nZ;
figure(2)
for t=1:1:nTime
    figure(4)
    errorbar(zShow,mu_H_uq1(zShow,t),var_H_uq1(zShow,t),'-')
    hold on
    plot(zShow,mid_H_uq1(zShow,t),'-')
    
    errorbar(zShow,mu_H_uq2(zShow,t),var_H_uq2(zShow,t),'--')
    plot(zShow,mid_H_uq2(zShow,t),'--')
    hold off
    title(sprintf('Mean Variance and Median @t=%i',t))
%     legend('All KL basis','Truncation KL basis')
    drawnow
%     frame(t)=getframe;    %comment to save cpu time and memory
end










end
