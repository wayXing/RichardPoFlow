% function permeaFieldApprox_Demo()
% Author:   Wei Xing
% History:  25/07/2017  file created

%% set parameter
lengthX=40;
deltaX=4;
nX=lengthX/deltaX+1;

lengthY=40;
deltaY=4;
nY=lengthY/deltaY+1;

nSample=10;
muY=10;
lengthcale=10;
DeviationRatio=0.05;
nKL=10;

%% Main
[X,Y] = ndgrid(0:deltaX:lengthX,0:deltaY:lengthY);
location=[X(:),Y(:)];


% K=permeaField(location,lengthcale,muY,DeviationRatio,nSample);
K=permeaFieldApprox(location,lengthcale,muY,DeviationRatio,nSample,nKL);


K=reshape(K,nX,nY,nSample);


for i=1:nSample
    figure(i)

    pcolor(X,Y,K(:,:,i))
    shading interp;
    colormap jet;
    colorbar
    title(sprintf('Permeability field'))
    
%     contourf(X,Y,K(:,:,i))
% %     colormap(hot)
%     shading interp;
%     colorbar
%     title(sprintf('Permeability field'))
    
end