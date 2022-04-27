% Validation from sythetic data
restoredefaultpath;clc;clear;close all
addpath(genpath([pwd '\functions']));
[Maps,alldata] = Calibration_2DKIII(3,1,5);
[J,KI,KII,KIII] = Abaqus_2D_KIII(Maps);

%% HR-EBSD example, you will need mtex and Strain2Disp_FE
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
[Maps.X1,Maps.Y1,Maps.Ux,Maps.Uy,Maps.M4] = FE_OOM([...
    alldata(:,1:2) alldata(:,4:5) alldata(:,7)],'Linear',Maps.results);
% 3D 2nd
percentTry = 100:-10:10;    ctx = 1;   % to reduce data density when data is huge
while ctx<length(percentTry)+1
    try
        [inData,~] = smotherData(alldata,percentTry(ctx));
        X = inData(:,1);	X(isnan(inData(:,4)))=NaN;  inData(:,1) = X;
        FE_OOM(inData,'Linear',Maps.results);
        ctx = length(percentTry)+1;
    catch err
        ctx = ctx+1;
        warning(err.message);
        if ctx == length(percentTry)+1
            rethrow(err);
        end
    end
end
[X1,Y1,Maps.Z1,~,~,Maps.Uz] = plo3dUxy([Maps.results '\3D_Integrated_Uxy']);
Uz = mean(Maps.Uz,3);
F11 = scatteredInterpolant(X1(:),Y1(:),Uz(:),'natural');
% Evaluate FE data on grid of experimental results
newM = F11(Maps.X1(:),Maps.Y1(:));
Maps.Uz = reshape(newM,size(Maps.Ux,1),size(Maps.Ux,2));
save(Maps.SavingD,'Maps','-append')
cd(PWD)
% Calc. SIFs
[J,KI,KII,KIII] = Abaqus_2D_KIII(Maps);