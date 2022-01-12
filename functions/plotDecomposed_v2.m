function plotDecomposed_v2(Maps)
Um3 = 1/2*(Maps.Uz-flipud(Maps.Uz));
        % in case it is zero as Abaqus won't work
Um4 = 1/2*(Maps.Uz+flipud(Maps.Uz)) + ones(size(Maps.Ux))*1e-12;
        
%%
s1=subplot(2,3,1);  	contourf(Maps.X1,Maps.Y1,Maps.Ux,'LineStyle','none'); 	
title('U_x','fontsize',20);
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(2,3,2);  	contourf(Maps.X1,Maps.Y1,Maps.Uy,'LineStyle','none'); 	
title('U_y','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(2,3,3);  	contourf(Maps.X1,Maps.Y1,Maps.Uz,'LineStyle','none'); 	
title('U_z','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;

s5=subplot(2,3,6);  	contourf(Maps.X1,Maps.Y1,Um3,'LineStyle','none'); 	
title('U^{III}_{Asy.}','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;     colorbar off;
s6=subplot(2,3,5);  	contourf(Maps.X1,Maps.Y1,Um4,'LineStyle','none'); 	
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off; 
title('U^{III}_{Sy.}','fontsize',20);
addScale([2 3 6],[Maps.X1(:) Maps.Y1(:)]);
%
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'southoutside','position', [0.1356 0.2432 0.2151 0.0322] );
h.Label.String = [ 'U [\mum]']; 
set([s1 s2 s3 s5 s6],"clim",caxis);
%}
set(gcf,'position',[1 51 1460 950]);  