restoredefaultpath;clc;clear;close all
addpath(genpath([pwd '\functions']));
[Maps,M4] = Calibration_2DKIII(3,1,5);
[~,J,KI,KII,KIII] = Abaqus_2D_KIII_v2(Maps,M4);
%{
X = 0;
for iO = 0.31:-Maps.Maps.stepsize:-0.31
    Y = 0; X = X+1;
    for iV = 0.31:-Maps.Maps.stepsize:-0.31
        Y = Y+1;
    Maps.xo = [iO;-0.99];	Maps.yo = [iV;iV];
    x(X,Y) = iO;        y(X,Y) = iV;
    try
        [~,J,KI,KII,KIII] = Abaqus_2D_KIII(Maps,M4);
        
        rmdir('Abaqus Output','s');
    catch
        J.true = NaN;	J.div = NaN;    J.K.true = NaN;	J.K.div = NaN;
        KI.true = NaN;  KI.div = NaN;   KII.true = NaN;	KII.div = NaN;
        KIII.true = NaN;KIII.div = NaN;
    end
    j(X,Y) = J.true;        je(X,Y) = J.div;
    jk(X,Y) = J.K.true;     jke(X,Y) = J.K.div;
    kI(X,Y) = KI.true;    	kIe(X,Y) = KI.div;
    kII(X,Y) = KII.true;  	kIIe(X,Y) = KII.div;
    kIII(X,Y) = KIII.true; 	kIIIe(X,Y) = KIII.div;
    clear J KI KII KIII
    end
end
save('A:\OneDrive - Nexus365\GitHub\Abaqus_2D_KIII\Error.mat');
%}
%{
%%
clc;clear;close all
set(0,'defaultAxesFontSize',20);       set(0,'DefaultLineMarkerSize',20)
load('A:\OneDrive - Nexus365\GitHub\Abaqus_2D_KIII\Error.mat',...
    'j','jk','kI','kII','kIII','je','jke','kIe','kIIe','kIIIe','x','y');
J = 1-diag(j/j(8,8));    Je = diag(je/j(8,8));
JK = 1-diag(jk/jk(8,8));    JKe = diag(jke/jk(8,8));
KI = 1-diag(kI/kI(8,8));    KIe = diag(kIe/kI(8,8));
KII = 1-diag(kII/kII(8,8));    KIIe = diag(kIIe/kII(8,8));
KIII = 1-diag(kIII/kIII(8,8));    KIIIe = diag(kIIIe/kIII(8,8));

set(0,'defaultAxesFontSize',23); 	set(0,'DefaultLineMarkerSize',20)
subplot(2,3,1);contourf(x-0.02,y-0.02,(jk-jk(8,8))/jk(8,8),'LineStyle','none');
c=colorbar; axis image; c.Label.String = 'J_{Err}';
hold on; line([0,0],[-0.3,0.29],'color','w','LineStyle','--','linewidth',2);
line([-0.3,0.29],[0,0],'color','w','LineStyle','--','linewidth',2);hold off
subplot(2,3,2);contourf(x-0.02,y-0.02,1-kI/kI(8,8),'LineStyle','none');
c=colorbar; axis image; c.Label.String = 'KI_{Err}';axis off
subplot(2,3,4);contourf(x-0.02,y-0.02,1-kII/kII(8,8),'LineStyle','none');
c=colorbar; axis image; c.Label.String = 'KII_{Err}';axis off
subplot(2,3,5);contourf(x-0.02,y-0.02,1-kIII/kIII(8,8),'LineStyle','none');
c=colorbar; axis image; c.Label.String = 'KIII_{Err}';axis off
hold on; line([0.3,-0.3],[-0.3,0.3],'color','w','LineStyle','--','linewidth',2); hold off
brighten(0.5);brighten(0.5);colormap jet
fig=subplot(1,3,3);%fig=figure;
set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
yyaxis left;    hold on;
errorbar(diag(x),KI(:),KIe(:),'--.k','LineWidth',2,'MarkerFaceColor','k','markersize',25);
errorbar(diag(x),KII(:),KIIe(:),'--.b','MarkerEdgeColor','b','LineWidth',2,'markersize',25);
errorbar(diag(x),KIII(:),KIIIe(:),'--.c','LineWidth',2,'MarkerFaceColor','k','markersize',25);
ylabel('K_{Err}'); hold off
 ylim([-1.5 1.5]);
yyaxis right;
errorbar(diag(x),JK(:),JKe(:),'--.r','LineWidth',2,'MarkerFaceColor','r','markersize',25);
hold on; line([0,0],[-30,30],'color','k','LineStyle',':','linewidth',2);
line([-0.4,0.4],[0,0],'color','k','LineStyle',':','linewidth',2);hold off
ylabel('J_{Err}');        ylim([-0.7 0.7]);
xlabel('Distance from the tip (mm)');
[~,objh]=legend('K_{I}','K_{II}','K_{III}','J_{int}','location','best');
objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 25);
xlim([-0.33 0.33]);
set(gcf,'position',[-78 60 2034 883]);
saveas(gcf, 'Error.fig');
saveas(gcf, 'Error.tif');    close
%}
%{
%%
clc;clear
[Maps,M4] = Calibration_2DKIII(3,1,5);
count=0;
noise.X=[0 0.0001 0.0004 0.001 0.004 0.01 0.04 0.1 0.4 1 4 10 40 100];
for iX=noise.X
    % Create the noise values that we'll add to a.
    M5 = M4;
    M5.Ux = M4.Ux+randn(size(M4.Ux)) .* iX * 1e-6;
    M5.Uy = M4.Uy+randn(size(M4.Ux)) .* iX * 1e-6;
    M5.Uz = M4.Uz+randn(size(M4.Ux)) .* iX * 1e-6;
    [~,J,KI,KII,KIII] = Abaqus_2D_KIII(Maps,M5); % as desigignated maps
    fclose all; rmdir('Abaqus Output','s');
    count = count+1;
    noise.j(count) = J.K.true;           noise.je(count) = J.K.div;
    noise.ki(count) = KI.true;          noise.kie(count) = KI.div;
    noise.kii(count) = KII.true;         noise.kiie(count) = KII.div;
    noise.kiii(count) = KIII.true;        noise.kiiie(count) = KIII.div;
end
save('A:\OneDrive - Nexus365\GitHub\Abaqus_2D_KIII\Error.mat','noise','-append');
%
%%
set(0,'defaultAxesFontSize',18);       set(0,'DefaultLineMarkerSize',15)
load('A:\OneDrive - Nexus365\GitHub\Abaqus_2D_KIII\Error.mat','noise');
noise.X(noise.X>1.1)=NaN;
close all;
f=@(a,x) a*x; u = -6:0.1:0;
subplot(2,2,1);errorbar(noise.X,1-noise.ki(:)/noise.ki(1),noise.kie(:)/noise.ki(1),'--.k','LineWidth',2,'markersize',25);
xlabel('Induced Noise (%)');ylabel('%K_{I}^{Err}');  grid on; xlim([0.00008 1.4])
set(gca, 'XScale', 'log');xticks([1e-5 1e-4 1e-3 1e-2 0.1 1]); ylim([-0.5 0.5])

subplot(2,2,2);errorbar(noise.X,1-noise.kii(:)/noise.kii(1),noise.kiie(:)/noise.kii(1),'--.k','markersize',25,'LineWidth',2);
xlabel('Induced Noise (%)');ylabel('%K_{II}^{Err}'); grid on; xlim([0.00008 1.4])
set(gca, 'XScale', 'log');xticks([1e-5 1e-4 1e-3 1e-2 0.1 1]); ylim([-1 1])

subplot(2,2,3);errorbar(noise.X,1-noise.kiii(:)/noise.kiii(1),noise.kiiie(:)/noise.kiii(1),'--.k','LineWidth',2,'markersize',25);
xlabel('Induced Noise (%)');ylabel('%K_{III}^{Err}');  grid on; xlim([0.00008 1.4])
set(gca, 'XScale', 'log');xticks([1e-5 1e-4 1e-3 1e-2 0.1 1]); ylim([-0.1 0.1])

subplot(2,2,4);errorbar(noise.X,1-noise.j(:)/noise.j(1),noise.je(:)/noise.j(1),'--.k','LineWidth',2,'markersize',25);
xlabel('Induced Noise (%)');ylabel('%J_{int}^{Err}');  grid on; xlim([0.00008 1.4])
set(gca, 'XScale', 'log');xticks([1e-5 1e-4 1e-3 1e-2 0.1 1]); ylim([-0.15 0.15])

set(gcf,'position',[169 111 1100 730]);
saveas(gcf,'noise.tif');
saveas(gcf,'noise.fig'); close
%}