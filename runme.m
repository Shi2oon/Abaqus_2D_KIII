clear all; close all; clc

%% Input
% Material properties.
E = 1e3;                                                                    % Young's Modulus
v = 0.3;                                                                    % Poisson ratio
k = (3-v)/(1+v);                                                            % Bulk modulus
G = E/(2*(1 + v));                                                          % Shear modulus
% SIF loading
KI = 1;                                                                     % Mode I SIF
KII = 1;                                                                    % Mode II SIF
KIII = 1;                                                                   % Mode III SIF                                       
% Domain size (square, crack tip at centre).
sz = 50;
% Noise
nszsmr=5;                                                                   % Smeering size
ncnout=[0.25,0.25,0];                                                   	% Outlier position
nszout=0.02;                                                                % Outlier size
namout=0.01;                                                             	% Outlier amplitude
nrand=2e-5;                                                                 % Random noise

%% Anayltical displacement data.
GridSpacing = 1/sz*2;
lin = GridSpacing*(ceil(-1/GridSpacing)+1/2):GridSpacing:GridSpacing*(floor(1/GridSpacing)-1/2);
DataSize = [numel(lin),numel(lin),1];
[pX,pY,pZ] = meshgrid(lin,lin,0);
[th,r] = cart2pol(pX,pY);
% Displacement data
uX = 1/2*KI/G*sqrt(r/(2*pi)).*(+cos(th/2).*(k-cos(th)))+...
	1/2*KII/G*sqrt(r/(2*pi)).*(+sin(th/2).*(k+2+cos(th)));
uY = 1/2*KI/G*sqrt(r/(2*pi)).*(+sin(th/2).*(k-cos(th)))+...
    1/2*KII/G*sqrt(r/(2*pi)).*(-cos(th/2).*(k-2+cos(th)));
uZ = 2*KIII/G*sqrt(r/(2*pi)).*sin(th/2);
% Mask
mask=ones(DataSize);
mask(floor(DataSize(1)/2),1:floor(DataSize(2)/2))=nan;
% Add noise
ROI = ones(DataSize); ROI([1,end],:)=0; ROI(:,[1,end])=0;
roiids = find(ROI==0);
% Smeering
nXsmr = smoothdata(uX,'gaussian',nszsmr)-uX; nXsmr(roiids) = 0;
nYsmr = smoothdata(uY,'gaussian',nszsmr)-uY; nYsmr(roiids) = 0;
nZsmr = smoothdata(uZ,'gaussian',nszsmr)-uZ; nZsmr(roiids) = 0;
% Outlier
g = exp(-((pX-ncnout(1)).^2+(pY-ncnout(2)).^2+(pZ-ncnout(3)).^2)./(2*nszout^2));
nXout = g*namout; nXout(roiids) = 0;
nYout = g*namout; nYout(roiids) = 0;
nZout = g*namout; nZout(roiids) = 0; 
% Random noise
nXrnd = randn(DataSize)*nrand; nXrnd(roiids) = 0;
nYrnd = randn(DataSize)*nrand; nYrnd(roiids) = 0;
nZrnd = randn(DataSize)*nrand; nZrnd(roiids) = 0;
% Noisy dispalcement data
uXn = uX+round(nXsmr+nXout+nXrnd,5);
uYn = uY+round(nYsmr+nYout+nYrnd,5);
uZn = uZ+round(nZsmr+nZout+nZrnd,5);

%% Outlier removal
%uXm=uX; uYm=uY; uZm=uZ;
[uXm,uYm,uZm] = removeoutliers('grubbs',0.01,mask,uXn,uYn,uZn,'Tol',0.01);
%[uXm,uYm,uZm] = removeoutliers('median',3,mask,uXn,uYn,uZn,'Tol',0.01);

%% Show data
cntrs=linspace(min([uX(:);uY(:);uZ(:)]),max([uX(:);uY(:);uZ(:)]),20);
figure;
subplot(1,3,1); contour(pX,pY,uX.*mask,cntrs,'k-'); hold on; axis equal;
subplot(1,3,1); contour(pX,pY,uXn.*mask,cntrs,'r-'); 
subplot(1,3,1); contour(pX,pY,uXm.*mask,cntrs,'b-'); 
ylabel('y'); xlabel('x'); title('u_X')
subplot(1,3,2); contour(pX,pY,uY.*mask,cntrs,'k-'); hold on; axis equal;
subplot(1,3,2); contour(pX,pY,uYn.*mask,cntrs,'r-'); 
subplot(1,3,2); contour(pX,pY,uYm.*mask,cntrs,'b-');
xlabel('x'); title('u_Y')
subplot(1,3,3); contour(pX,pY,uZ.*mask,cntrs,'k-'); hold on; axis equal;
subplot(1,3,3); contour(pX,pY,uZn.*mask,cntrs,'r-'); 
subplot(1,3,3); contour(pX,pY,uZm.*mask,cntrs,'b-');
xlabel('x'); ; title('u_Z')
legend('Analytical','Noisy','Outlier removal')
%% SIF extraction - Decomposition method.
% Mode I–III Decomposition of the J-integral from DIC Displacement Data, Strain(2015), 51, 492-503
% Mode I
uXd(:,:,1) = 1/2*(uXm+flipud(uXm));
uYd(:,:,1) = 1/2*(uYm-flipud(uYm));
uZd(:,:,1) = zeros(DataSize);
% Mode II
uXd(:,:,2) = 1/2*(uXm-flipud(uXm));
uYd(:,:,2) = 1/2*(uYm+flipud(uYm));
uZd(:,:,2) = zeros(DataSize);
% Mode III
uXd(:,:,3) = zeros(DataSize);
uYd(:,:,3) = zeros(DataSize);
uZd(:,:,3) = 1/2*(uZm-flipud(uZm));
% Loop through decomposed displacement data
for m = 1:3
    % Strain
    [uXXd(:,:,m),uXYd(:,:,m),uXZd(:,:,m)] = crackgradient(uXd(:,:,m),GridSpacing);
    [uYXd(:,:,m),uYYd(:,:,m),uYZd(:,:,m)] = crackgradient(uYd(:,:,m),GridSpacing);
    [uZXd(:,:,m),uZYd(:,:,m),uZZd(:,:,m)] = crackgradient(uZd(:,:,m),GridSpacing);
    % Infitisimal strain tensor
    eXXd = uXXd;
    eYYd = uYYd;
    eZZd = uZZd;
    eXYd = 1/2*(uXYd+uYXd);
    eXZd = 1/2*(uXZd+uZXd);
    eYZd = 1/2*(uYZd+uZYd);
    % Chauchy stress tensor, assming linear-elastic, isotropic material
    sXXd = E/(1-v^2)*(eXXd+v*(eYYd+eZZd));
    sYYd = E/(1-v^2)*(eYYd+v*(eXXd+eZZd));
    sZZd = E/(1-v^2)*(eZZd+v*(eXXd+eYYd));
    sXYd = 2*G*eXYd;
    sXZd = 2*G*eXZd;
    sYZd = 2*G*eYZd;
    % Elastic strain energy
    Wd = 0.5*(eXXd.*sXXd+eYYd.*sYYd+eZZd.*sZZd+2*eXYd.*sXYd+2*eXZd.*sXZd+2*eYZd.*sYZd);
end
% Generate q field
celw = 5;  % Width of area contour. Has to be an odd number.
dQdX = ones(DataSize)/(GridSpacing*celw);
dQdX = flipud(tril(flipud(tril(dQdX))))-flipud(triu(flipud(triu(dQdX))));
dQdY = dQdX';
dQdX(dQdX.*dQdY~=0) = 0;
% Domain integral
dA = ones(DataSize).*GridSpacing^2;
JAd = ((sXXd.*uXXd+sXYd.*uYXd+sXZd.*uZXd-Wd).*dQdX+(sYYd.*uYXd+sXYd.*uXXd+sYZd.*uZXd).*dQdY).*dA;
JAd(isnan(JAd))=0;
% Contour selection
mid = floor(DataSize(1)/2);
[a,b] = meshgrid(1:DataSize(1));
linecon=round(max((abs(a-mid-1/2)),abs(b-mid-1/2)));
% Summation of integrals
for ii = 1:mid-1
    areaID = linecon>=(ii-floor(celw/2)) & linecon<=(ii+floor(celw/2)); 
    Jd(:,ii,:) = sum(sum(JAd.*areaID));
end
% Equivalent SIF
Kd(1:2,:) = sqrt(Jd(1:2,:)*E);	% Mode I & II
Kd(3,:) = sqrt(Jd(3,:)*2*G);	% Mode III
% Plot contours
figure;
plot(Kd','-o'); hold on; grid on; axis tight;
plot([ones(3,1),ones(3,1)*size(Kd,2)]',repmat([KI;KII;KIII],1,2)','-k');
legend('KI','KII','KIII')
ylabel('SIF'); xlabel('Contour (#)');
 
%% Support function
function [cx,cy,cz]=crackgradient(c,dx)
c=squeeze(c);
[row,~]=size(c);
midr=floor(row/2);
ctop=c(1:midr,:);
cbot=c(midr+1:end,:);
[cxtop,cytop]=gradient(ctop,dx);
[cxbot,cybot]=gradient(cbot,dx);
cx=[cxtop;cxbot];
cy=[cytop;cybot];
cz=zeros(size(cx));
end