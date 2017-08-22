% Porous network generation 2D Demo script
%
% Equation:
%
% Input parameters:
%
% Output parameters:
%
% Examples: see Demo
%
% % Short description of example, followed by Matlab code line
% >> matlab code of example
%
% See also: 
% Author:   Wei Xing
% History:  21/08/2017  file created
%
clear
%% Setup
Dim=2;              %Domain dimensionality
nNode=1000;

lengthX=20;
lengthY=20;
length=[lengthX,lengthY]'; 

lengthScaleX=100;     %use for determining the effective distance. set larger number for higher correlation (less randomness)
lengthScaleY=1;       %use for determining the effective distance. set larger number for higher correlation (less randomness)
lengthcale=[lengthScaleX,lengthScaleY]'; 

% maxEdgeDistance=1;  
maxEdgeDistance=sqrt(lengthX.^2+lengthY.^2)/(nNode^(1/Dim));     %maximum allowed edge distance. adjust according to lengthX/Y and nNode


muY=1;                           %mean of pore body
deviationRatio=0.5;
sigmaY=(muY*deviationRatio)^2;   %variance of pore body 

seed=101;

%% Generate node location
iMethod=2;
switch iMethod  
    case 1  %Grid locations     
        [X,Y] = ndgrid(0:1:lengthX,0:1:lengthY);
        X=[X(:),Y(:)];
        nNode=size(X,1);
        
    case 2 % random locations
        X = lhsdesign(nNode,Dim);
%         X(:,1)=X(:,1).*lengthX;
%         X(:,2)=X(:,2).*lengthY;
        X=X*diag(length);
        
end
% scatter(X(:,1),X(:,2));


%% Generate pore body

%Define covariance function (also means the variogram)
iMethod=2;
switch iMethod  
    case 1      % Isotropic Gaussian(exp) covariance structure
        distance = pdist(X);
        distanceMatrix = squareform(distance);
        H=distanceMatrix./lengthcale;       %H is the effective seperation distance (Matrix)
        SigmaY=exp(-H) .*sigmaY;
          
    case 2      % anisotropic Gaussian(exp) covariance structure
%         lengthcale=[lengthScaleX,lengthScaleY]'; 
%         distance = pdist(X);
        eDistance = pdist(X*spdiags(1./lengthcale,0,Dim,Dim)); % effective seperation distance
        H=squareform(eDistance);    %H is the effective seperation distance Matrix form
        SigmaY=exp(-H) .*sigmaY;
   
end


[muX,SigmaX]=LogN2N(muY*ones(nNode,1),SigmaY);

nSample=1;
% rng(seed);  %pseudo random
poreBody = exp( mvnrnd( muX , SigmaX , nSample ))';

% scatter(X(:,1),X(:,2),poreBody.*100);


%% Network connection using Delaunay Triangulation
DT = delaunayTriangulation(X);

triplot(DT)
% triplot(DT,X(:,1),X(:,2))
hold on
scatter(X(:,1),X(:,2),poreBody.*200);
hold off

%Trim the connections

%get edge list
edgeList = edges(DT);
%get edge length
distance = squareform (pdist(X));
edgeLinearInd = sub2ind([nNode,nNode],edgeList(:,1),edgeList(:,2));
edgeDist = distance(edgeLinearInd);
%Remove high distance edge
edgeList=edgeList(edgeDist<=maxEdgeDistance,:);
%form adjacency matrix
Adj = sparse([edgeList(:, 1),edgeList(:, 2)], [edgeList(:, 2),edgeList(:, 1)],1,nNode,nNode);

%% Plot
figure(2)
% gplot (Adj, X, 'o-',  'MarkerFaceColor','g','MarkerSize',poreBody.*100);
gplot (Adj, X, '+-');
hold on 
scatter(X(:,1),X(:,2),poreBody.*200);
% scatter(X(:,1),X(:,2),(poreBody.*5).^2);
hold off

G = graph(Adj);
% G.Nodes=X;
figure(3)
p=plot(G)
p.XData=X(:,1);
p.YData=X(:,2);



