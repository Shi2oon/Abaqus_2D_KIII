function Va = ConverXCourt2xEBSD(DirxEBSD)
Maps=load(DirxEBSD);
% rotation
Va.W11 = Maps.rmap_w11;	    Va.W12 = Maps.rmap_w12;      Va.W13 = Maps.rmap_w13;
Va.W21 = Maps.rmap_w21;	    Va.W22 = Maps.rmap_w22;      Va.W23 = Maps.rmap_w23;
Va.W31 = Maps.rmap_w31;	    Va.W32 = -Maps.rmap_w23;      Va.W33 = Maps.rmap_w33;
Va.PH = Maps.rmap_geompkhgt;      Va.MAE = Maps.rmap_mae;

% stress
Va.S11 = Maps.rmap_s11;	    Va.S12 = Maps.rmap_s12;      Va.S13 = Maps.rmap_s31;
Va.S21 = Maps.rmap_s12;	    Va.S22 = Maps.rmap_s22;      Va.S23 = Maps.rmap_s23;
Va.S31 = Maps.rmap_s31;	    Va.S32 = Maps.rmap_s23;      Va.S33 = Maps.rmap_s33;
%strain
Va.E11 = Maps.rmap_e11;	    Va.E12 = Maps.rmap_e12;      Va.E13 = Maps.rmap_e31;
Va.E21 = Maps.rmap_e12;	    Va.E22 = Maps.rmap_e22;      Va.E23 = Maps.rmap_e23;
Va.E31 = Maps.rmap_e31;	    Va.E32 = Maps.rmap_e23;      Va.E33 = Maps.rmap_e33;
% displacement gradient tensor
Va.A11 = Va.E11+Va.W11; Va.A12 = Va.E12+Va.W12; Va.A13 = Va.E13+Va.W13;
Va.A21 = Va.E21+Va.W21; Va.A22 = Va.E22+Va.W22; Va.A23 = Va.E23+Va.W23;
Va.A31 = Va.E31+Va.W31; Va.A32 = Va.E32+Va.W32; Va.A33 = Va.E33+Va.W33;
%
Va.Stiffness  = Maps.stiffnessvalues;
Va.X   = Maps.xpos;        Va.Y   = Maps.ypos;
Va.Version = 'xCourt';
Va.RefID = Maps.grain_number;
Va.GrainData = NaN(size(Va.E11));

% stepsize
uko = unique(Va.X );                Va.stepsize =(abs(uko(1)-uko(2)));
Va.A11 = Va.E11+Va.W11; Va.A12 = Va.E12+Va.W12; Va.A13 = Va.E13+Va.W13;
Va.A21 = Va.E21+Va.W21; Va.A22 = Va.E22+Va.W22; Va.A23 = Va.E23+Va.W23;
Va.A31 = Va.E31+Va.W31; Va.A32 = Va.E32+Va.W32; Va.A33 = Va.E33+Va.W33;
Va.Wo = (1/2).*(Va.S11.*Va.E11 + Va.S12.*Va.E12 + Va.S13.*Va.E13 +...
    Va.S21.*Va.E21 + Va.S22.*Va.E22 + Va.S23.*Va.E23 +...
    Va.S31.*Va.E31 + Va.S32.*Va.E32 + Va.S33.*Va.E33);
Va.nu  =  Va.Stiffness(1,2)/(Va.Stiffness(1,1)+ Va.Stiffness(1,2));
Va.E   =  Va.Stiffness(1,1)*(1-2*Va.nu)*(1+Va.nu)/(1-Va.nu);
Va.units.xy = 'um';       Va.units.S  = 'GPa';      Va.units.W = 'rad';
Va.units.E  = 'Abs.';     Va.units.St = 'GPa';

try
    load([DirxEBSD '_GND'],'Total_GND_Density_',...
        'Total_Edge_Dislocation_Density_','Total_Screw_Dislocation_Density_');
    Va.GND = Total_GND_Density_;
    Va.Edge = Total_Edge_Dislocation_Density_;
    Va.Screw = Total_Screw_Dislocation_Density_;
catch err
    disp(err.message)
    disp('GND is randomised')
    Va.GND = randi(1e15,size(Va.A21));
end

%%
if ~exist([erase(DirxEBSD,'.ctf') '_EBSD.mat'],'file')
    clc;    warning('off'); tic; close all
    cd('C:\Users\ak13\OneDrive - National Physical Laboratory\GitHub\MTEX_Pipeline')
    addpath([pwd '\mtex-5.2.beta2'])
    % addpath('P:\Abdo\GitHub\MyEBSD\Routine');   startup_mtex;
    set(0, 'DefaultFigureVisible', 'on');       set(0,'defaultAxesFontSize',25)
    addpath([pwd '\Routine'])
    path = erase(DirxEBSD, '.ctf');         mkdir(path);
    ebsd = loadEBSD([erase(DirxEBSD, '.ctf') '.ctf'],'interface','ctf');%...
    %                                 ,'convertEuler2SpatialReferenceFrame');


    CS = ebsd.CSList;
    for iv = 1:length(ebsd.indexedPhasesId)
        if  strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Ferrite, bcc (New)') ||...
                strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Iron-alpha')
            CS{ebsd.indexedPhasesId(iv)}.opt.type='bcc';
            ebsd.CSList{ebsd.indexedPhasesId(iv)}.mineral = 'Ferrite';
            CS{ebsd.indexedPhasesId(iv)}.mineral = 'Ferrite';
        elseif  strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Austenite, fcc (New)') ||...
                strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Iron-Austenite')
            CS{ebsd.indexedPhasesId(iv)}.opt.type='fcc';
            ebsd.CSList{ebsd.indexedPhasesId(iv)}.mineral = 'Austenite';
            CS{ebsd.indexedPhasesId(iv)}.mineral = 'Austenite';
        elseif  strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Silicon')
            CS{ebsd.indexedPhasesId(iv)}.opt.type='fcc';
        elseif      strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Moissanite 6H')
            CS{ebsd.indexedPhasesId(iv)}.opt.type='fcc';
            ebsd.CSList{ebsd.indexedPhasesId(iv)}.mineral = 'Moissanite';
            CS{ebsd.indexedPhasesId(iv)}.mineral = 'Moissanite';
        elseif  strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Moissanite')
            CS{ebsd.indexedPhasesId(iv)}.opt.type='fcc';
        elseif  strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Ni-superalloy') ||...
                strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Ni') ||...
                strcmpi(CS{ebsd.indexedPhasesId(iv)}.mineral,'Nickel-cubic')
            CS{ebsd.indexedPhasesId(iv)}.opt.type='fcc';
            ebsd.CSList{ebsd.indexedPhasesId(iv)}.mineral = 'Nickel-cubic';
            CS{ebsd.indexedPhasesId(iv)}.mineral = 'Nickel-cubic';
        end
    end
    setMTEXpref('xAxisDirection','east');% west for bruker
    %         setMTEXpref('yAxisDirection','south');
    setMTEXpref('zAxisDirection','intoPlane');% outOfPlane for bruker
    [ebsd,CS,grains,~,Misanglede] = GBs(ebsd,path,CS);
    [SF,sS,SFM,mP,TraceSy] = SchmidTaylorEBSD(CS,ebsd,Misanglede,grains,path);
    save([erase(DirxEBSD,'.ctf') '_EBSD.mat'],'ebsd','CS','grains',...
        'SF','sS','SFM','mP','TraceSy');
else
    setMTEXpref('xAxisDirection','east');% west for bruker
    %         setMTEXpref('yAxisDirection','south');
    setMTEXpref('zAxisDirection','intoPlane');% outOfPlane for bruker
    load([erase(DirxEBSD,'.ctf') '_EBSD.mat'],'ebsd','CS','grains')

end
[~,~,ebsd_line] = LineEBSD(ebsd);
[C]=StiffnessEBSD(ebsd_line{2});
Va.Stiffness = C.Korsunsky;
Va.R = C.R;
Va.Mat = ebsd_line{2}.mineral;
Va.E = C.E;
Va.nu = C.nu;
Va.G = C.G;

%%
    function [lineData,lengthslope,ebsd_line]=LineEBSD(ebsd)
        % plot
        close all
        [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',5*degree);
        ebsd(grains(grains.grainSize<=20))   = [];
        [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',5*degree);
        %                     grains          = grains('indexed');
        %                     grains.boundary = grains.boundary('indexed');
        %                     grains          = smooth(grains,5);
        hold on;    plot(grains.boundary);     hold off

        for i=1:length(ebsd.indexedPhasesId)
            plot(ebsd(ebsd.mineralList{ebsd.indexedPhasesId(i)}),...
                ebsd(ebsd.mineralList{ebsd.indexedPhasesId(i)}).orientations); hold on
            gB = grains.boundary(ebsd.mineralList{ebsd.indexedPhasesId(i)}...
                ,ebsd.mineralList{ebsd.indexedPhasesId(i)});
            [xi,~]=size(gB);
            if xi~=0
                plot(gB,gB.misorientation.angle./degree,'linewidth',1);
            end
        end
        set(gcf,'position',[500,100,950,700]); hold off; colormap(jet(256));

        % grap
        title('select a point inside the grain you are targeting?')
        [xdata, ydata] = ginput(2); hold on;
        [StepSize] = CalcStepSize(ebsd);
        % Vertical line segment

        lineSec =  [xdata(1)   ydata(1); xdata(2) ydata(2)];%x2=x1 && y2>y1
        ebsd_line{1} = spatialProfile(ebsd,lineSec);       % line esbd data
        line(lineSec(:,1),lineSec(:,2),'linewidth',2,'LineStyle','--','Color','k');

        lineSec =  [xdata(1)+StepSize   ydata(1); xdata(2)+StepSize ydata(2)];%x2=x1 && y2>y1
        ebsd_line{2} = spatialProfile(ebsd,lineSec);       % line esbd data
        line(lineSec(:,1),lineSec(:,2),'linewidth',1,'LineStyle',':','Color','k');

        lineSec =  [xdata(1)-StepSize   ydata(1); xdata(2)-StepSize ydata(2)];%x2=x1 && y2>y1
        ebsd_line{3} = spatialProfile(ebsd,lineSec);       % line esbd data
        line(lineSec(:,1),lineSec(:,2),'linewidth',1,'LineStyle',':','Color','k'); hold off

        % find location in the data set
        xplot=ebsd_line{1}.prop.x;     yplot=ebsd_line{1}.prop.y;
        lineData=ebsd_line{1}.rotations.angle./degree;

        for i=1:length(xplot)
            lengthslope(i)=sqrt(abs(yplot(i)-yplot(1))^2+abs(xplot(i)-xplot(1))^2);
        end

    end


%% StiffnessEBSD
    function [C]=StiffnessEBSD(ebsd)
        % a code to calculate stifness matrix for each grain
        % you will need MTEX
        % there os an example included .. feel free to run it and to compare the
        % results to 'What is the Young’s Modulus of Silicon?', equation 8
        % https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5430873
        % ref. http://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_with_solutions.pdf
        CS = ebsd.CSList{ebsd.phase+1};
        csGlaucophane = crystalSymmetry(CS.pointGroup,'mineral',CS.mineral);
        % getting the stifness values, make sure the name is consistent with the
        % function
        stiffvalues   = lStiffness(CS.mineral);%in GPa
        % define the tensor
        C.R = ebsd.orientations.matrix;
        C.R = mean(C.R,3);
        % other options xEBSDV3, Dunne_StiffnessRot
        C.Korsunsky = Korsunsky_StiffnessRot(C.R,stiffvalues);%in GPa

        CGlaucophane = stiffnessTensor(C.Korsunsky,csGlaucophane);
        C.E  = CGlaucophane.YoungsModulus(vector3d.X);
        C.nu = CGlaucophane.PoissonRatio(vector3d.X,vector3d.Z); % in plane
        C.G  = CGlaucophane.shearModulus(vector3d.X,vector3d.Y); %in-plane shear
    end

    function DA=Korsunsky_StiffnessRot(R,L)
        %% Salvati, E., Sui, T. and Korsunsky, A.M., 2016.
        % Uncertainty quantification of residual stress evaluation by the FIB-DIC
        % ring-core method due to elastic anisotropy effects.
        % International Journal of Solids and Structures, 87, pp.61-69.
        R11 = R(1,1);   R12 = R(1,2);   R13 = R(1,3);
        R21 = R(2,1);   R22 = R(2,2);   R23 = R(2,3);
        R31 = R(3,1);   R32 = R(3,2);   R33 = R(3,3);
        TsA = [   R11^2,   R21^2,   R31^2,         2*R21*R31,         2*R11*R31,         2*R11*R21
            R12^2,   R22^2,   R32^2,         2*R22*R32,         2*R12*R32,         2*R12*R22
            R13^2,   R23^2,   R33^2,         2*R23*R33,         2*R13*R33,         2*R13*R23
            R12*R13, R22*R23, R32*R33, R22*R33 + R23*R32, R12*R33 + R13*R32, R12*R23 + R13*R22
            R11*R13, R21*R23, R31*R33, R21*R33 + R23*R31, R11*R33 + R13*R31, R11*R23 + R13*R21
            R11*R12, R21*R22, R31*R32, R21*R32 + R22*R31, R11*R32 + R12*R31, R11*R22 + R12*R21];

        TeA = [R11^2,     R21^2,     R31^2,           R21*R31,           R11*R31,           R11*R21
            R12^2,     R22^2,     R32^2,           R22*R32,           R12*R32,           R12*R22
            R13^2,     R23^2,     R33^2,           R23*R33,           R13*R33,           R13*R23
            2*R12*R13, 2*R22*R23, 2*R32*R33, R22*R33 + R23*R32, R12*R33 + R13*R32, R12*R23 + R13*R22
            2*R11*R13, 2*R21*R23, 2*R31*R33, R21*R33 + R23*R31, R11*R33 + R13*R31, R11*R23 + R13*R21
            2*R11*R12, 2*R21*R22, 2*R31*R32, R21*R32 + R22*R31, R11*R32 + R12*R31, R11*R22 + R12*R21];
        DA=inv(TsA)*L*TeA;
    end

    function [stiffvalues] = lStiffness(materials)
        %New code for stiffness matrices - CMM 2019. Gives the correct stiffness matrix for the phases. Most data either from xebsdv2,
        %or from crosscourt materials.txt (probably Simons and Wagg).
        phasestiffness(1).name='Silver';
        phasestiffness(1).stiff=[123.99,93.67,	93.67,	0,  0,	0;
            93.67,	123.99,	93.67,  0,	0,	0;
            93.67,	93.67, 123.99,	0,	0,	0;
            0,	0,	0,46.12, 0,	0;
            0,	0,  0,	0,46.12,  0;
            0,  0,	0,	0,	0,46.12];

        phasestiffness(2).name='Corundum';
        phasestiffness(2).stiff=[497.350000000000,163.970000000000,112.200000000000,-23.5800000000000,0,0;
            163.970000000000,497.350000000000,112.200000000000,23.5800000000000,0,0;
            112.200000000000,112.200000000000,499.110000000000,0,0,0;
            -23.5800000000000,23.5800000000000,0,147.390000000000,0,0;
            0, 0,0,0,147.390000000000,-23.5800000000000;
            0, 0,0,0,-23.5800000000000,166.690000000000];

        phasestiffness(3).name='Aluminium';
        phasestiffness(3).stiff=    [106.750000000000,60.4100000000000,60.4100000000000,0,0,0;
            60.4100000000000,106.750000000000,60.4100000000000,0,0,0;
            60.4100000000000,60.4100000000000,106.750000000000,0,0,0;
            0,0,0,28.3400000000000,0,0;
            0,0,0,0,28.3400000000000,0;
            0,0,0,0,0,28.3400000000000];

        phasestiffness(4).name='Diamond';
        phasestiffness(4).stiff=[1040,170,170,0,0,0;
            170,1040,170,0,0,0;
            170,170,1040,0,0,0;
            0,0,0,550,0,0;
            0,0,0,0,550,0;
            0,0,0,0,0,550];

        phasestiffness(5).name='Chromium';
        phasestiffness(5).stiff=[339.800000000000,58.6000000000000,58.6000000000000,0,0,0;
            58.6000000000000,339.800000000000,58.6000000000000,0,0,0;
            58.6000000000000,58.6000000000000,339.800000000000,0,0,0;
            0,0,0,99,0,0;
            0,0,0,0,99,0;
            0,0,0,0,0,99];

        phasestiffness(6).name='Copper';
        phasestiffness(6).stiff=[168.300000000000,122.100000000000,122.100000000000,0,0,0;
            122.100000000000,168.300000000000,122.100000000000,0,0,0;
            122.100000000000,122.100000000000,168.300000000000,0,0,0;
            0,0,0,75.7000000000000,0,0;
            0,0,0,0,75.7000000000000,0;
            0,0,0,0,0,75.7];

        phasestiffness(7).name='Lithium';
        phasestiffness(7).stiff=[13.5,	11.44,	11.44,	0	,0	,0	;
            11.44	,13.5	,11.44	,0	,0	,0	;
            11.44	,11.44	,13.5	,0	,0	,0	;
            0	,0	,0	,8.78	,0	,0	;
            0	,0	,0	,0	,8.78	,0	;
            0	,0	,0	,0	,0	,8.78];

        phasestiffness(8).name='Iron-alpha';
        phasestiffness(8).stiff=[230,   135,    135,    0,  0,  0;
            135,   230,    135,    0,  0,  0;
            135,   135,    230,    0,  0,  0;
            0,     0,      0,      117,0,  0;
            0,     0,      0,      0,  117,0;
            0,     0,      0,      0,  0,  117];

        phasestiffness(9).name='Iron-Austenite';
        phasestiffness(9).stiff=[231.4, 134.7,  134.7,  0,      0,      0;
            134.7, 231.4,  134.7,  0,      0,      0;
            134.7, 134.7,  231.4,  0,      0,      0;
            0,     0,      0,      116.4,  0,      0;
            0,     0,      0,      0,      116.4,  0;
            0,     0,      0,      0,      0,      116.4];

        phasestiffness(10).name='Iron-beta';
        %no stiffness matrix
        phasestiffness(11).name='Iron-delta';
        %no stiffness matrix

        phasestiffness(12).name='Hematite';
        phasestiffness(12).stiff=[242.430000000000,54.6400000000000,15.4200000000000,-12.4700000000000,0,0;
            54.6400000000000,242.430000000000,15.4200000000000,12.4700000000000,0,0;
            15.4200000000000,15.4200000000000,227.340000000000,0,0,0;
            -12.4700000000000,12.4700000000000,0,85.6900000000000,0,0;
            0,0,0,0,85.6900000000000,-12.4700000000000;
            0,0,0,0,-12.4700000000000,93.8950000000000];

        phasestiffness(13).name='Magnetite';
        phasestiffness(13).stiff=[273,106,106,0,0,0;
            106,273,106,0,0,0;
            106,106,273,0,0,0;
            0,0,0,97.1000000000000,0,0;
            0,0,0,0,97.1000000000000,0;
            0,0,0,0,0,97.1000000000000];

        phasestiffness(14).name='Magnesium';
        phasestiffness(14).stiff=[59.5000000000000,26.1200000000000,21.8000000000000,0,0,0;
            26.1200000000000,59.5000000000000,21.8000000000000,0,0,0;
            21.8000000000000,21.8000000000000,61.5500000000000,0,0,0;
            0,0,0,16.3500000000000,0,0;
            0,0,0,0,16.3500000000000,0;
            0,0,0,0,0,16.6900000000000];

        phasestiffness(15).name='Manganese-gamma';
        %no stiffness matrix

        phasestiffness(16).name='Molybdnenum-bcc';
        phasestiffness(16).stiff=[463.700000000000,157.800000000000,157.800000000000,0,0,0;
            157.800000000000,463.700000000000,157.800000000000,0,0,0;
            157.800000000000,157.800000000000,463.700000000000,0,0,0;
            0,0,0,109.200000000000,0,0;
            0,0,0,0,109.200000000000,0;
            0,0,0,0,0,109.200000000000];

        phasestiffness(17).name='Halite';
        phasestiffness(17).stiff=[49.4700000000000,12.8800000000000,12.8800000000000,0,0,0;
            12.8800000000000,49.4700000000000,12.8800000000000,0,0,0;
            12.8800000000000,12.8800000000000,49.4700000000000,0,0,0;
            0,0,0,12.8700000000000,0,0;
            0,0,0,0,12.8700000000000,0;
            0,0,0,0,0,12.8700000000000];

        phasestiffness(18).name='GaAs';
        phasestiffness(18).stiff=[48.8000000000000,41.4000000000000,41.4000000000000,0,0,0;
            41.4000000000000,48.8000000000000,41.4000000000000,0,0,0;
            41.4000000000000,41.4000000000000,48.8000000000000,0,0,0;
            0,0,0,14.8000000000000,0,0;
            0,0,0,0,14.8000000000000,0;
            0,0,0,0,0,14.8000000000000]; %DO CHECK WHETHER THESE ARE CORRECT - will depend on manufacturing route.

        phasestiffness(19).name='Nickel-cubic';
        phasestiffness(19).stiff=[249 152 152 0 0 0;
            152 249 152 0 0 0;
            152 152 249 0 0 0;
            0 0 0 124 0 0;
            0 0 0 0 124 0;
            0 0 0 0 0 124];

        phasestiffness(20).name='Olivine';
        phasestiffness(20).stiff=[324,59,79,0,0,0;
            59,198,78,0,0,0;
            79,78,249,0,0,0;
            0,0,0,66.7000000000000,0,0;
            0,0,0,0,81,0;
            0,0,0,0,0,79.3000000000000];

        phasestiffness(21).name='Quartz';
        phasestiffness(21).stiff=[86.8000000000000,7.04000000000000,11.9100000000000,-18.0400000000000,0,0;
            7.04000000000000,86.8000000000000,11.9100000000000,18.0400000000000,0,0;
            11.9100000000000,11.9100000000000,105.750000000000,0,0,0;
            -18.0400000000000,18.0400000000000,0,58.2000000000000,0,0;
            0,0,0,0,58.2000000000000,-18.0400000000000;
            0,0,0,0,-18.0400000000000,39.8800000000000];

        phasestiffness(22).name='Silicon'; %the orthotropic stifness matrix
        % ref: What is the Young's Modulus for Silicon?
        phasestiffness(22).stiff=  [165.7,  63.9,   63.9,   0,      0,      0;
            63.9,   165.7,  63.9,   0,      0,      0;
            63.9,   63.9,   165.7,  0,      0,      0;
            0,      0,      0,      79.6,   0,      0;
            0,      0,      0,      0,      79.6,   0;
            0,      0,      0,      0,      0,      79.6];

        phasestiffness(23).name='Titanium-alpha';
        phasestiffness(23).stiff=[162.4,92,69,0,0,0;
            92,162.4,69,0,0,0;
            69,69,180.7,0,0,0;
            0,0,0,46.7,0,0;
            0,0,0,0,46.7,0;
            0,0,0,0,0,35.2];


        phasestiffness(24).name='Titanium-beta'; %https://doi.org/10.3390/met8100822
        phasestiffness(24).stiff=[134,110,110,0,0,0;
            110,134,110,0,0,0;
            110,110,110,0,0,0;
            0,0,0,55,0,0;
            0,0,0,0,55,0;
            0,0,0,0,0,55];

        phasestiffness(25).name='Baddeleyite'; %WILLI PABST, GABRIELA TICHÁ, EVA GREGOROVÁ, 2004
        phasestiffness(25).stiff=[327,100,62,0,0,0;
            100,327,62,0,0,0;
            62,62,264,0,0,0;
            0,0,0,64,0,0;
            0,0,0,0,64,0;
            0,0,0,0,0,64];

        phasestiffness(26).name='Zirconium-bcc'; %"Lattice dynamics of hcp and bcc zirconium", Jerel L. Zarestky, 1979
        phasestiffness(26).stiff=[143.4,72.8,65.3,0,0,0;
            72.8,143.4,65.3,0,0,0;
            65.3,65.3,164.8,0,0,0;
            0,0,0,32,0,0;
            0,0,0,0,32,0;
            0,0,0,0,0,75.3];

        phasestiffness(27).name='Zirconium-hcp'; %xebsd2
        phasestiffness(27).stiff=[143.4,72.8,65.3,0,0,0;
            72.8,143.4,65.3,0,0,0;
            65.3,65.3,164.8,0,0,0;
            0,0,0,32,0,0;
            0,0,0,0,32,0;
            0,0,0,0,0,35.3];

        phasestiffness(28).name='Moissanite';
        phasestiffness(28).stiff=[504,  98,     56, 0	,0	,0	;
            98,   504,    56, 0	,0	,0	;
            56,   56,     566,0	,0	,0	;
            0,    0,      0,  170	,0	,0	;
            0,    0,      0,  0	,170,0	;
            0,    0,      0,  0	,0	,154];

        phasestiffness(29).name='Tungsten'; %crosscourt, not xebsd2
        phasestiffness(29).stiff=[522.39, 204.37, 204.37,   0,      0,      0;
            204.37, 522.39, 204.37,   0,      0,      0;
            204.37, 204.37, 522.39,   0,      0,      0;
            0	, 0,        0,      160.83, 0,      0;
            0	, 0,        0,      0,      160.83, 0;
            0	, 0,        0,      0,      0,      160.83];

        % Titanium aluminide (1/1)
        phasestiffness(30).name  = 'TiAl-Gamma';
        phasestiffness(30).stiff = [200	62	84	0	0	0
            62	200	84	0	0	0
            84	84	174	0	0	0
            0	0	0	112	0	0
            0	0	0	0	112	0
            0	0	0	0	0	41];

        % Aluminium titanium (1/3)
        phasestiffness(31).name  = 'Ti3Al';
        phasestiffness(31).stiff = [192	85	63	0	0	0
            85	192	63	0	0	0
            63	63	234	0	0	0
            0	0	0	60	0	0
            0	0	0	0	60	0
            0	0	0	0	0	54];

        % titanium Aluminium  (1/3)
        phasestiffness(32).name  = 'TiAl3';
        phasestiffness(32).stiff = [189	64	64	0	0	0
            64	189	64	0	0	0
            64	64	189	0	0	0
            0	0	0	73	0	0
            0	0	0	0	73	0
            0	0	0	0	0	73];

        % Modify names
        if  strcmpi(materials,'Ferrite, bcc (New)') || strcmpi(materials,'Ferrite')
            materials = 'Iron-alpha';
        elseif strcmpi(materials,'Austenite, fcc (New)') || strcmpi(materials, 'Austenite')
            materials = 'Iron-Austenite';
        elseif strcmpi(materials,'Moissanite 6H')
            materials = 'Moissanite';
        elseif strcmpi(materials,'Ni-superalloy')
            materials = 'Nickel-cubic';
        end

        %
        t = strcmp({phasestiffness(:).name}, materials); % find the field with the right name
        stiffvalues=phasestiffness(t).stiff; %use the field for the stiffness matrix

    end


%% calulcate effective E
    function [E,v,G,Co] = effectiveE_nu(C)
        % this function caclulates effective Youn modulus and Possion ratio for a
        % ansitropic material based on this paper
        % Reference: https://doi.org/10.3390/cryst8080307
        BV = (C(1,1)+2*C(1,2))/3;               % Voigt bulk modulus
        GV = (C(1,1)-C(1,2)+3*C(4,4))/5;        % Voigt shear modulus

        S = C^-1;
        BR = 1/(3*S(1,1)+6*S(1,2));             % Reuss bulk modulus
        GR = 5/(4*S(1,1)-4*S(1,2)+3*S(4,4));    % Reuss shear modulus

        B = (BR+BV)/2;                          % Hill’s average bulk modulus
        G = (GR+GV)/2;                          % Hill’s average shear modulus
        E = 9*B*G/(3*B+G);                      % Young’s modulus (E)
        v = (3*B-E)/(6*B);                      % Poisson’s ratio
        Co = [];
    end

end

