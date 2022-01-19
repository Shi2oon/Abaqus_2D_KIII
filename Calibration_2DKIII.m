function [Maps,alldata] = Calibration_2DKIII(KI,KII,KIII)
%% Input
              close all                   
% Domain size (square, crack tip at centre).
Maps.Mat          = 'Calibration';
Maps.type         = 'A'; % 'A' if u want to use anistropic matrix or 'E' for linear elastic
Maps.input_unit   = 'um';        % meter (m) or milmeter (mm) or micrometer(um);
Maps.units.xy     = Maps.input_unit; 
Maps.units.S      = 'Pa';      
Maps.units.St     = 'Pa'; 
Maps.pixel_size   = 1;           % if DIC values are in pixel, 1 if in physical units;
Maps.Dim          = '3D';        % handles  2D and 3D data with option
Maps.Operation    = 'xED';       % Strain, xED = xEBSD, DIC = Displacement
Maps.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
Maps.unique       = 'Calibration';
sz = 50;

switch Maps.input_unit
    case 'm'
        saf = 1;
    case 'mm'
        saf = 1e3;
    case 'um'
        saf = 1e6;
    case 'nm'
        saf = 1e9;
end

Maps.E = 210e9;             % Young's Modulus
Maps.nu = 0.3;             % Poisson ratio
G = Maps.E/(2*(1 + Maps.nu));  % Shear modulus
    switch Maps.stressstat
        case 'plane_strain'
            kappa = 3 - (4 .* Maps.nu); % [/]
        case 'plane_stress'
            kappa = (3 - Maps.nu)./(1 + Maps.nu); % [/]
    end                                                           % Bulk modulus                                                  
% SIF loading
KI = KI*1e6;                                                                     % Mode I SIF
KII = KII*1e6;                                                                    % Mode II SIF
KIII = KIII*1e6; 
%
Maps.Stiffness = [1/Maps.E          -Maps.nu/Maps.E     -Maps.nu/Maps.E 0 0 0
                 -Maps.nu/Maps.E        1/Maps.E        -Maps.nu/Maps.E 0 0 0
                 -Maps.nu/Maps.E    -Maps.nu/Maps.E         1/Maps.E    0 0 0
                  0 0 0     2*(1+Maps.nu)/Maps.E                          0 0
                  0 0 0 0       2*(1+Maps.nu)/Maps.E                        0
                  0 0 0 0 0         2*(1+Maps.nu)/Maps.E];
Maps.Stiffness = Maps.Stiffness^-1;
%}
Maps.SavingD = [pwd];
Maps.results = [pwd];

%% Anayltical displacement data.
Maps.stepsize = 1/sz*2;
lin = Maps.stepsize*(ceil(-1/Maps.stepsize)+1/2):Maps.stepsize:Maps.stepsize*(floor(1/Maps.stepsize)-1/2);
[Maps.X,Maps.Y,Maps.Z] = meshgrid(lin,lin,0);
[th,r] = cart2pol(Maps.X,Maps.Y);
DataSize = [numel(lin),numel(lin),1];
Maps.X1 = Maps.X*saf;
Maps.Y1 = Maps.Y*saf;
Maps.Z1 = Maps.Z*saf;
Maps.X = Maps.X*saf;
Maps.Y = Maps.Y*saf;
Maps.Z = Maps.Z*saf;
% displacement data
Maps.Ux = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+cos(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa+2+cos(th))))*saf;
Maps.Uy = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(-cos(th/2).*(kappa-2+cos(th))))*saf;
Maps.Uz = ( 2*KIII/G*sqrt(r/(2*pi)).*sin(th/2))*saf;

Maps.stepsize = Maps.stepsize*saf;
 
Maps.xo = [-0.01;-0.99]*saf;        Maps.xm = [0.01;-0.99]*saf;
Maps.yo = [0.0026;0.0026]*saf;      Maps.ym = [0.03;-0.03]*saf;

%% JMAN approach (without FEM) - Standard J-integral.
[Maps.E11,Maps.E12,Maps.E13] = crackgradient(Maps.Ux,Maps.stepsize);
[Maps.E21,Maps.E22,Maps.E23] = crackgradient(Maps.Uy,Maps.stepsize);
[Maps.E31,Maps.E32,Maps.E33] = crackgradient(Maps.Uz,Maps.stepsize);

%%
eXX = Maps.E11;
eYY = Maps.E22;
eZZ = Maps.E33;
eXY = 1/2*(Maps.E12+Maps.E21);
eXZ = 1/2*(Maps.E13+Maps.E31);
eYZ = 1/2*(Maps.E23+Maps.E32);

% Chauchy stress tensor, assuming linear-elastic, isotropic material
Maps.S11 = Maps.E/(1-Maps.nu^2)*(eXX+Maps.nu*(eYY+eZZ));
Maps.S22 = Maps.E/(1-Maps.nu^2)*(eYY+Maps.nu*(eXX+eZZ));
Maps.S33 = Maps.E/(1-Maps.nu^2)*(eZZ+Maps.nu*(eXX+eYY));
Maps.S12 = 2*G*eXY;
Maps.S13 = 2*G*eXZ;
Maps.S23 = 2*G*eYZ;
Maps.E12 = eXY;     Maps.E21 = eXY;
Maps.E13 = eXZ;     Maps.E31 = eXZ;
Maps.E32 = eYZ;     Maps.E23 = eYZ;
alldata = [Maps.X(:) Maps.Y(:) Maps.Z(:) Maps.E11(:) Maps.E22(:) Maps.E33(:)...
           Maps.E12(:) Maps.E13(:) Maps.E23(:)]; 
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
