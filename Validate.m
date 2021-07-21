clc;clear; restoredefaultpath;close all
[Maps,M4] = Calibration_2DKIII(3,1,5);
[~,J,KI,KII,KIII] = Abaqus_2D_KIII(Maps,M4);