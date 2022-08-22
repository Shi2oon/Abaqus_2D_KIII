Data = importdata(fullfile(fIno,f{iV},fstruct(iO).name));
            [~,Maps ] = reshapeData([Data.data(:,1:2) ...
                zeros(length(Data.data(:,1)),1) Data.data(:,4:6)]);
            Maps.X = Maps.X1;     Maps.Y = Maps.Y1;     Maps.Z = Maps.Z1;
            Maps.Ux = Maps.Ux;    Maps.Uy = Maps.Uy;    Maps.Uz = Maps.Uz;
            clear Data
            Maps.results = erase(fullfile(fIno,f{iV},fstruct(iO).name),'.dat');
            Maps.stepsize = Maps.Y1(2,2);
            Maps.type       = 'E';
            Maps.stressstat = 'plane_stress';
            Maps.Operation  = 'DIC';
            Maps.nu = 0.36;         
            Maps.units.xy = 'mm';
            Maps.E = 3e9;           
            Maps.Mat = 'PMMA';
            Maps.pixel_size=1;
            [J,KI,KII,KIII] = Abaqus_2D_KIII(Maps);