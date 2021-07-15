function [Maps,alldata] = Calibration_2DKIII(KI,KII,KIII)
%% Input
              close all                   
% Domain size (square, crack tip at centre).
Maps.Mat          = 'Calibration';
Maps.type         = 'A';
Maps.input_unit   = 'm';        % meter (m) or milmeter (mm) or micrometer(um);
Maps.units.xy     = Maps.input_unit; 
Maps.units.S      = 'Pa';      
Maps.units.St      = 'Pa'; 
Maps.pixel_size   = 1;           % if DIC values are in pixel, 1 if in physical units;
Maps.Dim          = '3D';        % handles  2D and 3D data with option
Maps.Operation    = 'xED';       % Strain, xED = xEBSD, DIC = Displacement
Maps.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
Maps.unique       = 'Calibration';
sz = 100;

switch Maps.units.xy
    case 'm'
        saf = 1;
    case 'mm'
        saf = 1e3;
    case 'um'
        saf = 1e6;
    case 'nm'
        saf = 1e9;
end


Maps.E = 210e9;                                                                    % Young's Modulus
Maps.nu = 0.3;                                                                    % Poisson ratio
    switch Maps.stressstat
        case 'plane_strain'
            kappa = 3 - (4 .* Maps.nu); % [/]
        case 'plane_stress'
            kappa = (3 - Maps.nu)./(1 + Maps.nu); % [/]
    end                                                           % Bulk modulus
G = Maps.E/(2*(1 + Maps.nu));                                                          % Shear modulus
% SIF loading
KI = KI*1e6;                                                                     % Mode I SIF
KII = KII*1e6;                                                                    % Mode II SIF
KIII = KIII*1e6;     
% Maps.Stiffness = [1/Maps.E          -Maps.nu/Maps.E     -Maps.nu/Maps.E 0 0 0
%                  -Maps.nu/Maps.E        1/Maps.E        -Maps.nu/Maps.E 0 0 0
%                  -Maps.nu/Maps.E    -Maps.nu/Maps.E         1/Maps.E    0 0 0
%                   0 0 0     2*(1+Maps.nu)/Maps.E                          0 0
%                   0 0 0 0       2*(1+Maps.nu)/Maps.E                        0
%                   0 0 0 0 0         2*(1+Maps.nu)/Maps.E];
% Maps.Stiffness = Maps.Stiffness^-1;
Maps.SavingD = [pwd '\KI-II-III'];   mkdir(Maps.SavingD);

%% Anayltical displacement data.
Maps.stepsize = 1/sz*2;
lin = Maps.stepsize*(ceil(-1/Maps.stepsize)+1/2):Maps.stepsize:Maps.stepsize*(floor(1/Maps.stepsize)-1/2);
[Maps.X,Maps.Y,Maps.Z] = meshgrid(lin,lin,0);
[th,r] = cart2pol(Maps.X,Maps.Y);
DataSize = [numel(lin),numel(lin),1];
Maps.M4.X = Maps.X*saf;
Maps.M4.Y = Maps.Y*saf;
Maps.M4.Z = Maps.Z*saf;
% displacement data
Maps.M4.Ux = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+cos(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa+2+cos(th))))*saf;
Maps.M4.Uy = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(-cos(th/2).*(kappa-2+cos(th))))*saf;
Maps.M4.Uz = ( 2*KIII/G*sqrt(r/(2*pi)).*sin(th/2))*saf;

Maps.stepsize = Maps.stepsize*saf;

%{
%%
switch Maps.units.xy
    case 'm'
        saf = 1;
    case 'mm'
        saf = 1e-3;
    case 'um'
        saf = 1e-6;
    case 'nm'
        saf = 1e-9;
end
Maps.units.St = 'Pa';        Maps.units.xy = 'm';

%% Anayltical displacement data.
M4=Maps.M4;
tmp = sortrows([M4.X(:) M4.Y(:) M4.Z(:) M4.Ux(:) M4.Uy(:) M4.Uz(:)],[3,1,2]);
[~,dataum ] = reshapeData(tmp);
Maps.X = mean(dataum.X1,3)*saf;     Maps.Y = mean(dataum.Y1,3)*saf;
Maps.Z = mean(dataum.Z1,3)*saf;     Maps.Ux = mean(dataum.Ux,3)*saf;
Maps.Uy = mean(dataum.Uy,3)*saf;    Maps.Uz = mean(dataum.Uz,3)*saf;
Maps.stepsize = min([diff(unique(Maps.X(:)))' diff(unique(Maps.Y(:)))'])*saf;
%}
%%
figure; s1=subplot(1,4,1);  	contourf(Maps.M4.X,Maps.M4.Y,Maps.M4.Ux ,'LineStyle','none'); 	
title('U_x [\mum]','fontsize',20);
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(1,4,2);  	contourf(Maps.M4.X,Maps.M4.Y,Maps.M4.Uy,'LineStyle','none'); 	
title('U_y [\mum]','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(1,4,3);  	contourf(Maps.M4.X,Maps.M4.Y,Maps.M4.Uz,'LineStyle','none'); 	
title('U_z [\mum]','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
s4=subplot(1,4,4);contourf(Maps.M4.X,Maps.M4.Y,sqrt(Maps.M4.Ux.^2 + ...
    Maps.M4.Uy.^2+Maps.M4.Uz.^2),'LineStyle','none')
axis image; axis off; colormap jet;  box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;
addScale([1 4 4],[Maps.M4.X(:) Maps.M4.Y(:)]);title('U_{Mag} [\mum]','fontsize',20);

cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'southoutside','position', [0.3513 0.0993 0.3 0.03] );
h.Label.String = [ 'U [\mum]']; 
set([s1 s2 s3 s4],"clim",caxis);

set(gcf,'position',[1 41 1900 900]);  
% saveas(gcf,[fileparts(file) '\2D_Maps.E12z.tif'],'tiffn');     
% saveas(gcf,[fileparts(file) '\2D_Maps.E12z.fig']);  close 
 
Maps.xo = [0.0026;-1.1716]*saf;     Maps.xm = [0.0155,-1.0477]*saf;
Maps.yo = [0,-0.0026]*saf;          Maps.ym = [0.0181,-0.0181]*saf;

%% JMAN approach (without FEM) - Standard J-integral.
[Maps.E11,Maps.E12,Maps.E13] = crackgradient(Maps.M4.Ux,Maps.stepsize);
[Maps.E21,Maps.E22,Maps.E23] = crackgradient(Maps.M4.Uy,Maps.stepsize);
[Maps.E31,Maps.E32,Maps.E33] = crackgradient(Maps.M4.Uz,Maps.stepsize);
% Infitisimal strain tensor
eXX = Maps.E11;
eYY = Maps.E22;
eZZ = Maps.E33;
eXY = 1/2*(Maps.E12+Maps.E21);
eXZ = 1/2*(Maps.E13+Maps.E31);
eYZ = 1/2*(Maps.E23+Maps.E32);
alldata = [Maps.X(:) Maps.Y(:) Maps.Z(:) Maps.E11(:) Maps.E12(:) Maps.E13(:)...
    Maps.E21(:) Maps.E22(:) Maps.E23(:) Maps.E31(:) Maps.E32(:) Maps.E33(:)]; 
% Chauchy stress tensor, assuming linear-elastic, isotropic material
Maps.S11 = Maps.E/(1-Maps.nu^2)*(eXX+Maps.nu*(eYY+eZZ));
Maps.S22 = Maps.E/(1-Maps.nu^2)*(eYY+Maps.nu*(eXX+eZZ));
Maps.S33 = Maps.E/(1-Maps.nu^2)*(eZZ+Maps.nu*(eXX+eYY));
Maps.S12 = 2*G*eXY;
Maps.S13 = 2*G*eXZ;
Maps.S23 = 2*G*eYZ;
% Elastic strain energy
Maps.Wo = 0.5*(eXX.*Maps.S11+eYY.*Maps.S22+eZZ.*Maps.S33+2*eXY.*Maps.S12+...
    2*eXZ.*Maps.S13+2*eYZ.*Maps.S23);
Maps.GND = zeros(size(Maps.E11));

%{
%% Generate q field
celw = 1;  % Width of area contour. Has to be an odd number.
dQdX = ones(DataSize)/(Maps.stepsize*celw);
dQdX = flipud(tril(flipud(tril(dQdX))))-flipud(triu(flipud(triu(dQdX))));
dQdY = dQdX';
dQdX(dQdX.*dQdY~=0) = 0;
% Domain integral
dA = ones(DataSize).*Maps.stepsize^2;
JA = ((Maps.S11.*Maps.E11+Maps.S12.*Maps.E21+Maps.S13.*Maps.E31-Maps.Wo).*dQdX...
     +(Maps.S22.*Maps.E21+Maps.S12.*Maps.E11+Maps.S23.*Maps.E31).*dQdY).*dA;
% Contour selection
mid = floor(DataSize(1)/2);
[a,b] = meshgrid(1:DataSize(1));
linecon=round(max((abs(a-mid-1/2)),abs(b-mid-1/2)));
% Summation of integrals
for ii = 1:mid-1
    areaID = linecon>=(ii-floor(celw/2)) & linecon<=(ii+floor(celw/2)); 
    J.Raw(ii) = sum(sum(JA.*areaID));
end
% Equivalent energy release rate
Ja = 1/Maps.E*(KI^2+KII^2)+1/(2*G)*KIII^2;

%% Decomposition method.
% Mode I-III Decomposition of the J-integral from DIC Displacement Data, Strain(2015), 51, 492-503
% Mode I
uXd(:,:,1) = 1/2*(Maps.M4.Ux+flipud(Maps.M4.Ux));
uYd(:,:,1) = 1/2*(Maps.M4.Uy-flipud(Maps.M4.Uy));
uZd(:,:,1) = 1/2*(Maps.M4.Uz+flipud(Maps.M4.Uz)); % uz-flipud(uz) = 0
% Mode II
uXd(:,:,2) = 1/2*(Maps.M4.Ux-flipud(Maps.M4.Ux));
uYd(:,:,2) = 1/2*(Maps.M4.Uy+flipud(Maps.M4.Uy));
uZd(:,:,2) = zeros(DataSize);
% Mode III
uXd(:,:,3) = zeros(DataSize);
uYd(:,:,3) = zeros(DataSize);
uZd(:,:,3) = 1/2*(Maps.M4.Uz-flipud(Maps.M4.Uz));
% Loop through decomposed displacement data
for m = 1:3
    % Strain
    [Maps.E11d(:,:,m),Maps.E12d(:,:,m),Maps.E13d(:,:,m)] = crackgradient(uXd(:,:,m),Maps.stepsize);
    [Maps.E21d(:,:,m),Maps.E22d(:,:,m),Maps.E23d(:,:,m)] = crackgradient(uYd(:,:,m),Maps.stepsize);
    [Maps.E31d(:,:,m),Maps.E32d(:,:,m),Maps.E33d(:,:,m)] = crackgradient(uZd(:,:,m),Maps.stepsize);
    % Infitisimal strain tensor
    eXXd = Maps.E11d;
    eYYd = Maps.E22d;
    eZZd = Maps.E33d;
    eXYd = 1/2*(Maps.E12d+Maps.E21d);
    eXZd = 1/2*(Maps.E13d+Maps.E31d);
    eYZd = 1/2*(Maps.E23d+Maps.E32d);
    % Chauchy stress tensor, assming linear-elastic, isotropic material for
    % plane stress conditions
Maps.S11d = Maps.E/(1-Maps.nu^2)*(eXXd+Maps.nu*(eYYd+eZZd));
Maps.S22d = Maps.E/(1-Maps.nu^2)*(eYYd+Maps.nu*(eXXd+eZZd));
Maps.S33d = Maps.E/(1-Maps.nu^2)*(eZZd+Maps.nu*(eXXd+eYYd));
    Maps.S12d = 2*G*eXYd;
    Maps.S13d = 2*G*eXZd;
    Maps.S23d = 2*G*eYZd;
    % Elastic strain energy
    Wd = 0.5*(eXXd.*Maps.S11d+eYYd.*Maps.S22d+eZZd.*Maps.S33d+2*eXYd.*Maps.S12d+2*eXZd.*Maps.S13d+2*eYZd.*Maps.S23d);
end
% Generate q field
celw = 1;  % Width of area contour. Has to be an odd number.
dQdX = ones(DataSize)/(Maps.stepsize*celw);
dQdX = flipud(tril(flipud(tril(dQdX))))-flipud(triu(flipud(triu(dQdX))));
dQdY = dQdX';
dQdX(dQdX.*dQdY~=0) = 0;
% Domain integral
dA = ones(DataSize).*Maps.stepsize^2;
JAd = ((Maps.S11d.*Maps.E11d+Maps.S12d.*Maps.E21d+Maps.S13d.*Maps.E31d-Wd).*dQdX+(Maps.S22d.*Maps.E21d+Maps.S12d.*Maps.E11d+Maps.S23d.*Maps.E31d).*dQdY).*dA;
% Contour selection
mid = floor(DataSize(1)/2);
[a,b] = meshgrid(1:DataSize(1));
linecon=round(max((abs(a-mid-1/2)),abs(b-mid-1/2)));
% Summation of integrals
for ii = 1:mid-1
    areaID = linecon>=(ii-floor(celw/2)) & linecon<=(ii+floor(celw/2)); 
    J.K.Raw(:,ii,:) = sum(sum(JAd.*areaID));
end
% Equivalent SIF
Kd(1:2,:) = sqrt(J.K.Raw(1:2,:)*Maps.E*saf);	% Mode I & II
Kd(3,:) = sqrt(J.K.Raw(3,:)*2*G*saf);      % Mode III
J.K.Raw = sum(J.K.Raw);
clear KI KII KIII
KI.Raw   = Kd(1,:);
KII.Raw  = Kd(2,:);
KIII.Raw  = Kd(3,:);

%%
    J.Raw     = J.Raw*saf;             % in J/m^2
    J.K.Raw   = J.K.Raw.*saf;           % in J/m^2
    KI.Raw    = KI.Raw*1e-6;     % in MPa
    KII.Raw   = KII.Raw*1e-6;    % in MPa
    KIII.Raw   = KIII.Raw*1e-6;    % in MPa
    
%%
contrs   = length(J.Raw);        contrs = contrs - round(contrs*0.4);
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(contrs:end))))))+2;
if dic<1;       dic = 1;    end
J.true   = round(mean(rmoutliers(J.Raw(contrs:end))),dic);
J.div    = round(std(rmoutliers(J.Raw(contrs:end)),1),dic);
J.K.true = round(mean(rmoutliers(J.K.Raw(contrs:end))),dic);
J.K.div  = round(std(rmoutliers(J.K.Raw(contrs:end)),1),dic);
KI.true  = round(mean(rmoutliers( KI.Raw (contrs:end))),dic);
KI.div   = round(std(rmoutliers( KI.Raw (contrs:end)),1),dic);
KII.true = round(mean(rmoutliers( KII.Raw (contrs:end))),dic);
KII.div  = round(std(rmoutliers( KII.Raw (contrs:end)),1),dic);
KIII.true = round(mean(rmoutliers( KIII.Raw (contrs:end))),dic);
KIII.div  = round(std(rmoutliers( KIII.Raw (contrs:end)),1),dic);

%%
Contour = (1:length(J.Raw))*Maps.stepsize;
fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
yyaxis left;    hold on;
plot(Contour,Kd(1,:),'k--o','MarkerEdgeColor','k','LineWidth',4);
plot(Contour,Kd(2,:),'k--s','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');
plot(Contour,Kd(3,:),'k--<','MarkerEdgeColor','k','LineWidth',4');
ylabel('K (MPa m^{0.5})'); hold off
if min(Kd(:))>0;     ylim([0 max(Kd(:))+min(Kd(:))/3]);      end
yyaxis right;
plot(Contour,J.K.Raw,'r--o','MarkerEdgeColor','r','LineWidth',4);
hold on
plot(Contour,J.Raw,'r--s','MarkerEdgeColor','r','LineWidth',4);
plot(Contour,Ja,'r--<','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
hold off
ylabel('J [J/m^2]');        ylim([0 max(J.Raw)+min(J.K.Raw)/4]);
xlabel('DEFT [\mum]'); % two xlabel, Domain extension from the tip
xlabel('Contour Number');
legend(['K_{I} = '        num2str(KI.true)  ' ± ' num2str(KI.div)  ' MPa\surdm' ],...
    ['K_{II} = '       num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
    ['K_{III} = '      num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
    ['J_{integral} = ' num2str(J.true)   ' ± ' num2str(J.div)   ' J/m^2'],...
    ['J_{K} = ' num2str(J.K.true)   ' ± ' num2str(J.K.div)   ' J/m^2'],...
    ['J_{true} = ' num2str(Ja(1))   ' J/m^2'],...
    'location','northoutside','box','off');
grid on; xlim([0 max(Contour)+Maps.stepsize])
set(fig,'position',[60,10,850,1050]);   box off;  
saveas(fig, [fileparts(Maps.SavingD) '\KI-II-III\J_KI_II_III_calib.fig']);
saveas(fig, [fileparts(Maps.SavingD) '\KI-II-III\J_KI_II_III_calib.tif']);    close

saveas(gcf,[fileparts(Maps.SavingD) '\DIC.tif']);
saveas(gcf,[fileparts(Maps.SavingD)'\DIC.fig']);     close all
%}
close all
end
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