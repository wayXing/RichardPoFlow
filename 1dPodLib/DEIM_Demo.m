%% DEIM_Demo
clear
% Solve y(x,mu)=(1-x)cos(3*pi*mu*(x+1))*e^(-(1+x)*mu);

X=linspace(-1,1,100)';
Mu=linspace(1,3*pi,100)';

for i=1:length(Mu)
    Y(:,i)=(1-X).*cos(3*pi*Mu(i)*(X+1)).*exp(-(1+X)*Mu(i));   
end

[U,S,V]=svd(Y);

% Ur=U(:,1:20);
% [phi,Uu,P] = DEIM(Ur);
% D=Ur*inv(P'*Ur);

[phi,Uu,P] = DEIM(U);
Pr=P(:,1:10);
Ur=U(:,1:10);
Dr=Ur*inv(Pr'*Ur);


%% Example plot
% mu=2;
% y=(1-X).*cos(3*pi*mu*(X+1)).*exp(-(1+X)*mu);
% y_approx= D*(P'*y);
% 
% figure
% plot(X,y,X,y_approx);

%% MSSE
% Mu_test=linspace(1,3*pi,1000)';
% for i=1:length(Mu_test)
%     Y(:,i)=(1-X).*cos(3*pi*Mu_test(i)*(X+1)).*exp(-(1+X)*Mu_test(i));   
% end
% Y_approx= D*(P'*Y);
% 
% n=900;
% figure
% plot(X,Y(:,n),X,Y_approx(:,n));
% 
% SE=(Y-Y_approx).^2;
% MSSE=mean(SE(:));

%%

Mu_test=linspace(1,3*pi,1000)';
for i=1:length(Mu_test)
    Y(:,i)=(1-X).*cos(3*pi*Mu_test(i)*(X+1)).*exp(-(1+X)*Mu_test(i));   
end

Dim_deim=5:5:100;

for j =1:length(Dim_deim)
    Pr=P(:,1:Dim_deim(j));
    Ur=U(:,1:Dim_deim(j));
    Dr=Ur*inv(Pr'*Ur);
    Y_approx= Dr*(Pr'*Y);
    SE=(Y-Y_approx).^2;
    MSSE(j)=mean(SE(:));   
end
figure
plot(Dim_deim,log(MSSE),'*-')
title('MSSE VS Dimension of DEIM')








