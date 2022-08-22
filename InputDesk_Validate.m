% Validation from sythetic data
restoredefaultpath;clc;clear;close all
addpath(genpath([pwd '\functions']));
[Maps,~] = Calibration_2DKIII(3,1,-5);
[J,KI,KII,KIII] = Abaqus_2D_KIII(Maps);