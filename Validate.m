restoredefaultpath;clc;clear;close all
addpath(genpath([pwd '\functions']));
[Maps,M4] = Calibration_2DKIII(3,1,5);
[J,KI,KII,KIII] = Abaqus_2D_KIII_v2(Maps,M4);

%% HR-EBSD example, you will need mtex
% restoredefaultpath;clc;clear;close all
% run('A:\OneDrive - Nexus365\GitHub\mtex-5.2.beta2\install_mtex.m')
clc;clear;close all; PWD = pwd;
addpath(genpath([PWD '\functions']));
filename = [PWD '\Si_indent_21_XEBSD'];
% get HR_EBSD data, for include 3D in the title for 3D data
[Maps,alldata] = GetGrainData(filename,'Si_indent_21_XEBSD_3D');
Maps.results = fileparts(Maps.SavingD);
% download the Strain2Disp_FE at https://github.com/Shi2oon/Strain2Disp_FE
cd('A:\OneDrive - Nexus365\GitHub\Strain2Disp_FE')
% 2D first
[Maps.X1,Maps.Y1,Maps.Ux,Maps.Uy,Maps.M4] = ...
    FE_OOM([alldata(:,1:2) alldata(:,4:5) alldata(:,7)],'Linear',Maps.results);
% 3D 2nd
[~,~,~,~,~,Maps.Z1,Maps.Uz] = FE_OOM(alldata,'Linear',Maps.results);
plo3dUxy([Maps.results '\3D_Integrated_Uxy']);
cd(PWD)
% Calc. SIFs
[J,KI,KII,KIII] = Abaqus_2D_KIII_v2(Maps);