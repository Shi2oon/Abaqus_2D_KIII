function [X1,Y1,Z1,Ux,Uy,Uz] = plo3dUxy(file)
load(file); close all;
set(0,'defaultAxesFontSize',22);       set(0,'DefaultLineMarkerSize',14) 
tmp = sortrows(OutData,[3,1,2]);
[~,dataum ] = reshapeData(tmp);
X3D = dataum.X1;      Y3D = dataum.Y1;          Z3D = dataum.Z1;
Ux3D = dataum.Ux;     Uy3D = dataum.Uy;         Uz3D = dataum.Uz;

%%   
     Plot3D(sqrt(Ux3D.^2+Uy3D.^2+Uz3D.^2),X3D,Y3D,Z3D,'\mum','U_{mag}'); 
     axis tight; axis image
     saveas(gcf,[fileparts(file) '\3D_Umag.tif'],'tiffn');     
	 saveas(gcf,[fileparts(file) '\3D_Umag.fig']);  close
     
     X1 = X3D(:,:,1);     Y1 = Y3D(:,:,1);          Z1 = Z3D(:,:,1);
     Ux = Ux3D(:,:,1);   Uy = Uy3D(:,:,1);        Uz = Uz3D(:,:,1);
     
%% data some time the other way around      
    s4=subplot(1,1,1);contourf(X1,Y1,sqrt(Ux.^2+Uy.^2+Uz.^2),'LineStyle','none')
    axis image; colormap jet;  box off; %set(gca,'Xdir','reverse')
    c  =colorbar;	cU(4,:) = c.Limits;     colorbar off; axis off;title('U_{Mag}');
    addScale([1 1 1],[X1(:) Y1(:)]); xlabel('X'); ylabel('Y');
    set(gcf,'position',[1 41 800 800]); 
    set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
    
    reply = questdlg_timer(10,'Is the plot is on the right axis?','Axis','Y','N','Y');
    if reply=='N'
        X1 = Y3D(:,:,1)';         Y1 = X3D(:,:,1)';         Z1 = Z3D(:,:,1)';
        Ux = max(Ux3D,[],3)';     Uy = max(Uy3D,[],3)';     Uz = max(Uz3D,[],3)';
    end
%     save(file,'X1','Y1','Z1','Ux','Uy','Uz','X3D','Y3D','Z3D','Ux3D','Uy3D','Uz3D','-append');

%%
	[~,ii] = reshapeStrainData(alldata); %strain
    ii.Exx(ii.Exx==0) = NaN;    ii.Eyy(isnan(ii.Exx)) = NaN;    ii.Ezz(isnan(ii.Exx)) = NaN;
    ii.Exy(isnan(ii.Exx)) = NaN;    ii.Exz(isnan(ii.Exx)) = NaN;    ii.Eyz(isnan(ii.Exx)) = NaN;
%{
if ~strcmpi(M4.ScaleYN, 'Y')
    DispData = [X1(:) Y1(:) Ux(:) Uy(:) Uz(:)];
    newM(:,1:2) = [ii.X1(:) ii.Y1(:) ];
    Fx = scatteredInterpolant(DispData(:,1),DispData(:,2),DispData(:,3),'natural');
    newM(:,3) = Fx(ii.X1(:),ii.Y1(:)); 
    Fy = scatteredInterpolant(DispData(:,1),DispData(:,2),DispData(:,4),'natural');
    newM(:,4) = Fy(ii.X1(:),ii.Y1(:)); 
    Fz = scatteredInterpolant(DispData(:,1),DispData(:,2),DispData(:,5),'natural');
    newM(:,5) = Fz(ii.X1(:),ii.Y1(:));
    [~,dataum] = reshapeData(newM(:,1:4)); %strain
    X1 = dataum.X1;         Ux = dataum.Ux;         Uy = dataum.Uy;  
    X1(isnan(ii.Ux))=NaN;   Ux(isnan(ii.Ux))=NaN;   Uy(isnan(ii.Ux))=NaN;
    [~,dataum] = reshapeData([newM(:,1:2) newM(:,4:5)]); %strain
    Y1 = dataum.Y1;         Uz = dataum.Uy;
    Y1(isnan(ii.Ux))=NaN;   Uz(isnan(ii.Ux))=NaN;
    save(file,'X1','Y1','Z1','Ux','Uy','Uz','-append');
end
    %}
%%
close;s=surf(X1,Y1,sqrt(Ux.^2+Uy.^2+Uz.^2)); s.EdgeColor = 'none';
xlabel('X [\mum]');     ylabel('Y [\mum]'); zlabel('U_{Mag} [\mum]'); 
set(gcf,'position',[1 41 800 800]); 
set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
saveas(gcf,[fileparts(file) '\Surf_Uxyz.tif'],'tiffn');     
saveas(gcf,[fileparts(file) '\Surf_Uxyz.fig']);  close 
clear cU
close;s=surf(X1,Y1,Uz); s.EdgeColor = 'none';
xlabel('X [\mum]'); ylabel('Y [\mum]');zlabel('U_{z} [\mum]'); 
set(gcf,'position',[1 41 800 800]); 
set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
saveas(gcf,[fileparts(file) '\Surf_Uz.tif'],'tiffn');     
saveas(gcf,[fileparts(file) '\Surf_Uz.fig']);  close 

s1=subplot(1,4,1);  	contourf(X1,Y1,Ux,'LineStyle','none'); 	
title('U_x','fontsize',20);set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(1,4,2);  	contourf(X1,Y1,Uy,'LineStyle','none'); 	
title('U_y','fontsize',20);set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(1,4,3);  	contourf(X1,Y1,Uz,'LineStyle','none'); 	
title('U_z','fontsize',20);set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
s4=subplot(1,4,4);contourf(X1,Y1,sqrt(Ux.^2+Uy.^2+Uz.^2),'LineStyle','none')
axis image; axis off; colormap jet;  box off; %set(gca,'Ydir','reverse')
c  =colorbar;	     colorbar off;%cU(4,:) = c.Limits;
addScale([1 4 4],[X1(:) Y1(:)]);title('U_{Mag}','fontsize',20);
set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
%
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'southoutside','position', [0.3513 0.0993 0.3 0.03] );
h.Label.String = [ 'U [\mum]']; 
set([s1 s2 s3 s4],"clim",caxis);
%}
set(gcf,'position',[1 41 1900 900]);  
saveas(gcf,[fileparts(file) '\2D_Uxyz.tif'],'tiffn');     
saveas(gcf,[fileparts(file) '\2D_Uxyz.fig']);  close 
    
%%
s1=subplot(3,3,1);  	contourf(ii.X1,ii.Y1,ii.Exx,'LineStyle','none');
title(['' char(949)  '_{xx}'],'fontsize',19);
axis image; axis off; colormap jet; box off; set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(3,3,2);  	contourf(ii.X1,ii.Y1,ii.Exy,'LineStyle','none');title(['' char(949)  '_{xy}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(3,3,3);  	contourf(ii.X1,ii.Y1,ii.Exz,'LineStyle','none');title(['' char(949)  '_{xz}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
s5=subplot(3,3,5);  	contourf(ii.X1,ii.Y1,ii.Eyy,'LineStyle','none');title(['' char(949)  '_{yy}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;
s6=subplot(3,3,6);  	contourf(ii.X1,ii.Y1,ii.Eyz,'LineStyle','none');title(['' char(949)  '_{yz}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(5,:) = c.Limits;    colorbar off; 
s9=subplot(3,3,9);  	contourf(ii.X1,ii.Y1,ii.Ezz,'LineStyle','none');title(['' char(949)  '_{zz}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off;
addScale([3 3 9],[X1(:) Y1(:)]);
 
s4=subplot(3,3,4);  	contourf(X1,Y1,Ux,'LineStyle','none'); 	title('U_x','fontsize',19);
axis image; axis off;  box off; colormap jet;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cu(1,:) = c.Limits;     colorbar off; 
s7=subplot(3,3,7);  	contourf(X1,Y1,Uy,'LineStyle','none'); 	title('U_y','fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cu(2,:) = c.Limits;     colorbar off;
s8=subplot(3,3,8);  	contourf(X1,Y1,Uz,'LineStyle','none'); 	title('U_z','fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cu(3,:) = c.Limits;    colorbar off; 
%
cu(abs(cu)==1)=0;
set([s4 s7 s8],"clim",[min(cu(:)) max(cu(:))]); 
subplot(3,3,4); c=colorbar;  c.Label.String = [ 'U (\mum)']; 
c.Position = [00.1236 0.1133 0.0112 0.4611];
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'westoutside','position', [0.9011 0.1211 0.0121 0.7533] );
h.Label.String = ['' char(949)  '']; 
set([s1 s2 s3 s5 s6 s9],"clim",caxis); 
%}
set(gcf,'position',[1 -41 1900 1000]); 
saveas(gcf,[fileparts(file) '\3D_UE.tif'],'tiffn');     
saveas(gcf,[fileparts(file) '\3D_UE.fig']);  close 
%{
%%
s1=subplot(3,3,1);  	imagesc(unique(ii.X1),unique(ii.Y1),ii.Ux);title('' char(949)  '_{11}');
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(3,3,2);  	imagesc(unique(ii.X1),unique(ii.Y1),ij.Ux);title('' char(949)  '_{12}');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(3,3,3);  	imagesc(unique(ii.X1),unique(ii.Y1),ij.Uy);title('' char(949)  '_{13}');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
s5=subplot(3,3,5);  	imagesc(unique(ii.X1),unique(ii.Y1),ii.Uy);title('' char(949)  '_{22}');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;
s6=subplot(3,3,6);  	imagesc(unique(ii.X1),unique(ii.Y1),ij.Uz);title('' char(949)  '_{23}');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;    colorbar off; 
s9=subplot(3,3,9);  	imagesc(unique(ii.X1),unique(ii.Y1),ii.Uz);title('' char(949)  '_{33}');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off;
addScale([3 3 9],[X1(:) Y1(:)]);
 
s4=subplot(3,3,4);  	imagesc(unique(X1),unique(Y1),Ux); 	title('U_x');
axis image; axis off;  box off; colormap jet;
c  =colorbar;	cu(1,:) = c.Limits;     colorbar off; 
s7=subplot(3,3,7);  	imagesc(unique(X1),unique(Y1),Uy); 	title('U_y');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cu(2,:) = c.Limits;     colorbar off;
s8=subplot(3,3,8);  	imagesc(unique(X1),unique(Y1),Uz); 	title('U_z');
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cu(3,:) = c.Limits;    colorbar off; 
%
cu(abs(cu)==1)=0;
set([s4 s7 s8],"clim",[min(cu(:)) max(cu(:))]); 
subplot(3,3,4); c=colorbar;  c.Label.String = [ 'U [\mum]']; 
c.Position = [00.1236 0.1133 0.0112 0.4611];
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'westoutside','position', [0.9011 0.1211 0.0121 0.7533] );
h.Label.String = '' char(949)  ''; 
set([s1 s2 s3 s5 s6 s9],"clim",caxis); 

set(gcf,'position',[1 41 1900 900]); 
saveas(gcf,[fileparts(file) '\3D_UEm.tif'],'tiffn');     
saveas(gcf,[fileparts(file) '\3D_UEm.fig']);  close 
%}
%%
s1=subplot(3,3,1);  	contourf(ii.X1,ii.Y1,ii.Exx,'LineStyle','none');title(['' char(949)  '_{xx}'],'fontsize',19);
axis image; axis off; colormap jet; box off; set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(3,3,2);  	contourf(ii.X1,ii.Y1,ii.Exy,'LineStyle','none');title(['' char(949)  '_{xy}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(3,3,3);  	contourf(ii.X1,ii.Y1,ii.Exz,'LineStyle','none');title(['' char(949)  '_{xz}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
s5=subplot(3,3,5);  	contourf(ii.X1,ii.Y1,ii.Eyy,'LineStyle','none');title(['' char(949)  '_{yy}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;
s6=subplot(3,3,6);  	contourf(ii.X1,ii.Y1,ii.Eyz,'LineStyle','none');title(['' char(949)  '_{yz}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(5,:) = c.Limits;    colorbar off; 
s9=subplot(3,3,9);  	contourf(ii.X1,ii.Y1,ii.Ezz,'LineStyle','none');title(['' char(949)  '_{zz}'],'fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off;
addScale([3 3 9],[X1(:) Y1(:)]);
 
s4=subplot(3,3,4);  	contourf(X1,Y1,Ux,'LineStyle','none'); 	title('U_x','fontsize',19);
axis image; axis off;  box off; colormap jet;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cu(1,:) = c.Limits;     colorbar off; 
s7=subplot(3,3,7);  	contourf(X1,Y1,Uy,'LineStyle','none'); 	title('U_y','fontsize',19);
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cu(2,:) = c.Limits;     colorbar off;
s8=subplot(3,3,8);  	contourf(X1,Y1,sqrt(Ux.^2+Uy.^2),'LineStyle','none'); 	title('U_{Mag}');
axis image; axis off; colormap jet; box off;set(gca,'Xdir','reverse');set(gca,'Ydir','reverse');
c  =colorbar;	cu(3,:) = c.Limits;    colorbar off; 
%
cu(abs(cu)==1)=0;
set([s4 s7 s8],"clim",[min(cu(:)) max(cu(:))]); 
subplot(3,3,4); c=colorbar;  c.Label.String = [ 'U (\mum)']; 
c.Position = [00.1236 0.1133 0.0112 0.4611];
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'westoutside','position', [0.9011 0.1211 0.0121 0.7533] );
h.Label.String = ['' char(949)  '']; 
set([s1 s2 s3 s5 s6 s9],"clim",caxis); 
%}
set(gcf,'position',[1 -41 1900 1000]); 
saveas(gcf,[fileparts(file) '\2D_UE.tif'],'tiffn');     
saveas(gcf,[fileparts(file) '\2D_UE.fig']);  close 
end