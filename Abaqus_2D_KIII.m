function [UnitOffset,J,KI,KII,KIII] = Abaqus_2D_KIII(Dir,M4)
% variable M4 is for the displacement components
%%
addpath([pwd '\functions'])
tmp = sortrows([M4.X(:) M4.Y(:) M4.Z(:) M4.Ux(:) M4.Uy(:) M4.Uz(:)],[3,1,2]);
% try
%     [Dir.RadEulerAng,Dir.rotCentre,tmp]=toShoeMake(tmp);
%     tmp(tmp(:,4)==0 & tmp(:,5)==0 & tmp(:,6)==0,:)=[];
%     saveas(gcf, [Dir.results '\Rot_Removed.tif']);    close all
% catch err
%     disp(err.message)
% end
[~,dataum ] = reshapeData(tmp);
Dir.Maps.X1 = mean(dataum.X1,3);    Dir.Maps.Y1 = mean(dataum.Y1,3);
Dir.Maps.Z1 = mean(dataum.Z1,3);    Dir.Maps.Ux = mean(dataum.Ux,3);
Dir.Maps.Uy = mean(dataum.Uy,3);    Dir.Maps.Uz = mean(dataum.Uz,3);

if strcmpi(Dir.type, 'A')
    [Dir.E,Dir.nu,Dir.G,Dir.Co] = effectiveE_v(Dir.Stiffness); % in Pa
else
    Dir.G = Dir.E/(2*(1 + Dir.nu));
end

if strcmpi(Dir.stressstat, 'plane_strain')
    Dir.E = Dir.E/(1-Dir.nu^2);% for HR-EBSD plane strain conditions
    Dir.G = Dir.E/(2*(1 + Dir.nu));
end

names = {'KI','KII','KIII'};
for iO=1:3
    Dirxyz = Dir;
    Dirxyz.unique = names{iO};
    if iO == 1      % Mode I
        Dirxyz.Maps.Ux = 1/2*(Dir.Maps.Ux+flipud(Dir.Maps.Ux));
        Dirxyz.Maps.Uy = 1/2*(Dir.Maps.Uy-flipud(Dir.Maps.Uy));
    elseif iO == 2 % Mode II
        Dirxyz.Maps.Ux = 1/2*(Dir.Maps.Ux-flipud(Dir.Maps.Ux));
        Dirxyz.Maps.Uy = 1/2*(Dir.Maps.Uy+flipud(Dir.Maps.Uy));
    elseif iO == 3 % Mode III
        Dirxyz.Maps.Ux = 1/2*(Dir.Maps.Uz-flipud(Dir.Maps.Uz));
        % in case it is zero as Abaqus won't work
        Dirxyz.Maps.Uy = 1/2*(Dir.Maps.Uz+flipud(Dir.Maps.Uz)) ... 
                         + ones(size(Dir.Maps.Ux))*1e-12;
    end
    alldata = Dirxyz.Maps;
    [DATA,UnitOffset,Dirxyz, Dirxyz.msk,SaveD] = ...
        Locate_Crack(alldata,Dirxyz.input_unit,Dirxyz.results,Dirxyz);
    % prepare and run abaqus cae
    [Abaqus,~] = PrintRunCode(Dirxyz, ...
        Dirxyz.msk,SaveD,ceil(min(size(DATA.X1))*0.5-2),UnitOffset);
%     UnitOffset = 1e-6;
    if iO == 1      % Mode I
%         Abaqus = [Dir.results '\Abaqus Output\KI'];
        [Jd,~,KI] = PlotKorJ(Abaqus,Dir.E,UnitOffset,1);
        loT(iO) = length(KI.Raw);
    elseif iO == 2 % Mode II
%         Abaqus = [Dir.results '\Abaqus Output\KII'];
        [Jd,~,~,KII] = PlotKorJ(Abaqus,Dir.E,UnitOffset,1);
        loT(iO) = length(KII.Raw);
    elseif iO==3 % fix KIII to shear rather than modulus
%         Abaqus = [Dir.results '\Abaqus Output\KIII'];
        [Jd,~,addKI,KIII] = PlotKorJ(Abaqus,Dir.E,UnitOffset,1);
        KIII.Raw = KIII.Raw*2*Dir.G/Dir.E;
        Jd.K.Raw = (KIII.Raw.*1e6).^2/(2*Dir.G);
        loT(iO)  = length(KIII.Raw);
    end
    JKRaw(iO,1:length(Jd.K.Raw)) = Jd.K.Raw;
    JRaw(iO,1:length(Jd.Raw)) = Jd.Raw;
end
JKRaw(4,1:length(addKI.Raw)) = -addKI.Raw;
JRaw(4,1:length(Jd.Raw)) = (-addKI.Raw.*1e6).^2/Dir.E;
J.JKIII = JKRaw;
J.JIII = JRaw;

%% Cut to the same contourf convergence
J.Raw = sum(J.JIII(:,1:min(loT)));
J.K.Raw = sum(J.JKIII(:,1:min(loT)));
J.Raw = J.Raw(1:min(loT));      
J.K.Raw = J.K.Raw(1:min(loT));
KI.Raw = KI.Raw(1:min(loT));    
addKI.Raw = -addKI.Raw(1:min(loT));
KII.Raw = KII.Raw(1:min(loT));
KIII.Raw = KIII.Raw(1:min(loT));

% add addtional KI to KI
KI.Raw = KI.Raw + addKI.Raw;

%%
contrs   = length(J.Raw);        contrs = contrs - round(contrs*0.4);
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(contrs:end))))))+2;
if dic<1;       dic = 1;    end
J.true   = round(mean(rmoutliers(J.Raw(contrs:end))),dic);
J.div    = round(std(rmoutliers(J.Raw(contrs:end)),1),dic);
J.K.true   = round(mean(rmoutliers(J.K.Raw(contrs:end))),dic);
J.K.div    = round(std(rmoutliers(J.K.Raw(contrs:end)),1),dic);

KI.true  = round(mean(rmoutliers(KI.Raw(contrs:end))),dic);
KI.div   = round(std(rmoutliers(KI.Raw(contrs:end)),1),dic);
KII.true = round(mean(rmoutliers(KII.Raw(contrs:end))),dic);
KII.div  = round(std(rmoutliers(KII.Raw(contrs:end)),1),dic);
KIII.true = round(mean(rmoutliers(KIII.Raw(contrs:end))),dic);
KIII.div  = round(std(rmoutliers(KIII.Raw(contrs:end)),1),dic);
%
%%
plotJKIII(KI,KII,KIII,J,Dir.Maps.stepsize,Dir.input_unit)
saveas(gcf, [Dir.results '\J_KI_II_III_abaqus.fig']);
saveas(gcf, [Dir.results '\J_KI_II_III_abaqus.tif']);    close all

save([Dir.results '\Abaqus_2D_KIII.mat'],'Dir','J','KI','KII','KIII','M4');

plotDecomposed(M4)
saveas(gcf, [Dir.results '\U_Dec.fig']);
saveas(gcf, [Dir.results '\U_Dec.tif']);    close
%}
end
