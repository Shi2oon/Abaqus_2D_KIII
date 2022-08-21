function [J,KI,KII,KIII] = Abaqus_2D_KIII(Maps)
% variable M4 is for the displacement components
%%
warning on; addpath([pwd '\functions'])
if size(Maps.Uz,3) ~=1
    Maps.Uz = Maps.Uz(:,:,1);
end
if size(Maps.Ux,2)*size(Maps.Ux,2) ~= size(Maps.Uz,2)*size(Maps.Uz,2)
    Maps.Uz = Maps.Uz(2:end,2:end);
end

%%
if strcmpi(Maps.type, 'A')
    [Maps.E,Maps.nu,Maps.G,Maps.Co] = effectiveE_v(Maps.Stiffness); % in Pa
else
    Maps.G = Maps.E/(2*(1 + Maps.nu));
end

if strcmpi(Maps.stressstat, 'plane_strain')
    Maps.E = Maps.E/(1-Maps.nu^2);% for HR-EBSD plane strain conditions
    Maps.G = Maps.E/(2*(1 + Maps.nu));
end

names = {'KI','KIII'};
for iO=1:2
    Dirxyz = Maps;
    Dirxyz.unique = names{iO};
    if iO == 1      % Mode I/II
        Dirxyz.Ux = Maps.Ux;
        Dirxyz.Uy = Maps.Uy;
    elseif iO == 2 % Mode III
        Dirxyz.Ux = Maps.Uz;%1/2*(Maps.Uz-flipud(Maps.Uz));
        % in case it is zero as Abaqus won't work
        Dirxyz.Uy = ones(size(Maps.Ux))*1e-12;%+ 1/2*(Maps.Uz+flipud(Maps.Uz)) ...
        
    end
    alldata = Dirxyz;
    [DATA,UnitOffset,Dirxyz, Dirxyz.msk,SaveD] = ...
        Locate_Crack(alldata,Dirxyz.units.xy,Dirxyz.results,Dirxyz);
    % prepare and run abaqus cae
    [Abaqus,~] = PrintRunCode(Dirxyz, ...
        Dirxyz.msk,SaveD,ceil(min(size(DATA.X1))*0.5-2),UnitOffset);
    
    %     UnitOffset = 1e-6;
    if iO == 1      % Mode I
        %         Abaqus = [Dir.results '\Abaqus Output\KI'];
        [Jd,~,KI,KII,Direction] = ...
            PlotKorJ(Abaqus,Maps.E,UnitOffset,1);
        if ~isempty(Direction.Raw)
            fprint('\nRecommended J-integral direction is %.2f',Direction.true)
            Ans = input('Do you want to adjust (N/Y)?');
            if strcmpi(Ans,'Y')
                [Abaqus] = Adjust4Direction(Abaqus,Direction.true);
                [Jd,~,KI,KII,Direction.Raw] = ...
                    PlotKorJ(Abaqus,Maps.E,UnitOffset,1);
            end
        end
        loT(iO) = length(KI.Raw);
        %         Maps.xo = Dirxyz.xo;        Maps.yo = Dirxyz.yo;
        %         Maps.xm = Dirxyz.xm;        Maps.ym = Dirxyz.ym;
        %         Maps.Operation = 'xED';
    elseif iO==2 % fix KIII to shear rather than modulus
        %         Abaqus = [Dir.results '\Abaqus Output\KIII'];
        [Jd,~,addKI,KIII,Dir] = PlotKorJ(Abaqus,Maps.E,UnitOffset,1);
        % correct from in-plane to out-of-plane shear
        KIII.Raw = KIII.Raw*2*Maps.G/Maps.E;
        %         Jd.Raw   = (KIII.Raw.*1e6).^2/(2*Dir.G);
        %         Jd.K.Raw = (KIII.Raw.*1e6).^2/(2*Dir.G);
        loT(iO)  = length(KIII.Raw);
        if ~isempty(Dir.Raw)
        Direction.Raw(iO,1:length(Dir))=Dir.Raw;
        Direction.true(iO)=Dir.true;
        Direction.div(iO)=Dir.div;
        end
    end
    % J when calculating the SIF (more accurate)
    JKRaw(iO,1:length(Jd.K.Raw)) = Jd.K.Raw;
    JRaw(iO,1:length(Jd.Raw)) = Jd.Raw; % J from J analysis
end
% JKRaw(3,1:length(addKI.Raw)) = (addKI.Raw.*1e6).^2/Dir.E;
% JRaw(3,1:length(addKI.Raw)) = (addKI.Raw.*1e6).^2/Dir.E; % J(I)r
J.JKIII = JKRaw;
J.JIII = JRaw;

%% Cut to the same contour convergence (IoT value)
J.Raw = sum(J.JIII(:,1:min(loT)));
J.K.Raw = sum(J.JKIII(:,1:min(loT)));

J.Raw = J.Raw(1:min(loT));
J.K.Raw = J.K.Raw(1:min(loT));
KI.Raw = KI.Raw(1:min(loT));
addKI.Raw = addKI.Raw(1:min(loT));
KI.addKI.Raw = KI.Raw + addKI.Raw ;% add addtional KI to KI
KII.Raw = KII.Raw(1:min(loT));
KIII.Raw = KIII.Raw(1:min(loT));

%%
contrs   = length(J.Raw);        contrs = contrs - round(contrs*0.4);
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(contrs:end))))))+2;
if dic<2;       dic = 2;    end
J.true   = round(mean(rmoutliers(J.Raw(contrs:end))),dic);
J.div    = round(std(rmoutliers(J.Raw(contrs:end)),1),dic);
J.K.true   = round(mean(rmoutliers(J.K.Raw(contrs:end))),dic);
J.K.div    = round(std(rmoutliers(J.K.Raw(contrs:end)),1),dic);
% J.addJ.Raw  = JRaw(3,1:min(loT));
% J.addJ.true = round(mean(rmoutliers(J.addJ.Raw(contrs:end))),dic);
% J.addJ.div  = round(std(rmoutliers(J.addJ.Raw(contrs:end)),1),dic);

KI.true  = round(mean(rmoutliers(KI.Raw(contrs:end))),dic);
KI.div   = round(std(rmoutliers(KI.Raw(contrs:end)),1),dic);
KI.addKI.true  = round(mean(rmoutliers(KI.addKI.Raw(contrs:end))),dic);
KI.addKI.div   = round(std(rmoutliers(KI.addKI.Raw(contrs:end)),1),dic);

KII.true = round(mean(rmoutliers(KII.Raw(contrs:end))),dic);
KII.div  = round(std(rmoutliers(KII.Raw(contrs:end)),1),dic);
KIII.true = round(mean(rmoutliers(KIII.Raw(contrs:end))),dic);
KIII.div  = round(std(rmoutliers(KIII.Raw(contrs:end)),1),dic);

%
%%
plotJKIII(KI,KII,KIII,J,Maps.stepsize,Maps.units.xy)
saveas(gcf, [Maps.results '\J_KI_II_III_abaqus.fig']);
saveas(gcf, [Maps.results '\J_KI_II_III_abaqus.tif']);    %close all

save([Maps.results '\Abaqus_2D_KIII.mat'],'Maps','J','KI',...
    'KII','KIII','Direction');

figure; plotDecomposed_v2(Maps)
saveas(gcf, [Maps.results '\U_Dec.fig']);
saveas(gcf, [Maps.results '\U_Dec.tif']);    %close
%}
end
