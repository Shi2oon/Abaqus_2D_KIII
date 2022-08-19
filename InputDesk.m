
[Maps,~] = Calibration_2DKIII(5,-2,-5);
% Calc. SIFs
[J,KI,KII,KIII] = Abaqus_2D_KIII(Maps);

%%
clc;clear
inO = [pwd '\S1_1206um_20kV_10nA_150uS_XEBSD_3D_Full_map'];
load(fullfile(inO, 'Crop & Rot Data.mat'),'Maps')
load(fullfile(inO, 'Linear_Integrated_Uxy.mat'),'Ux','Uy','X1','Y1')
Maps.results = inO;
[~,B] = fileparts(Maps.SavingD);
Maps.SavingD = fullfile(inO,B);
Maps.Ux = Ux;   Maps.Uy = Uy;   Maps.X1 = X1;   Maps.Y1 = Y1;
[X1,Y1,Maps.Z1,~,~,Maps.Uz] = plo3dUxy([Maps.results '\3D_Integrated_Uxy']);
Uz = mean(Maps.Uz,3);
F11 = scatteredInterpolant(X1(:),Y1(:),Uz(:),'natural');

% Evaluate FE data on grid of experimental results
newM = F11(Maps.X1(:),Maps.Y1(:));
Maps.Uz = reshape(newM,size(Maps.Ux,1),size(Maps.Ux,2));
save(Maps.SavingD,'Maps','-append')
% Calc. SIFs
[J,KI,KII,KIII] = Abaqus_2D_KIII(Maps);