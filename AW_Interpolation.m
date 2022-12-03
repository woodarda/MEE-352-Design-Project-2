function [Temperature, Pressure,vf, hfg, cpf, muf, kf, Prf] = AW_Interpolation(T)

%AW_Interpolation Summary of this function goes here
%   Detailed explanation goes here


Properties=readtable("WaterProperties.txt");
Properties=renamevars(Properties,["Var1","Var2", "Var3","Var4","Var5","Var6",...
    "Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16"],["Temp",...
    "Pressure","vf","vg","hfg","cpf","cpg","muf","mug","kf","kg","Prf",...
    "Prg","sigf","betaf","Temp2"]);

I=find(Properties.Temp==T);

if isempty(I)== 1  
I=find(Properties.Temp>T,1);

end
num1=Properties.Temp(I-1)-T;
den1=Properties.Temp(I-1)-Properties.Temp(I);
den2=Properties{I-1,:}-Properties{I,:};

Properties=Properties{I-1,:}-num1*den2/den1;
convToBaseUnits=[1 1 10^-3 1 10^3 10^3 10^3 10^-6 10^-6 10^-3 10^-3 1 1 ...
    10^-3 10^-6 1];
Properties=Properties.*convToBaseUnits;

Temperature=Properties(1);
Pressure=Properties(2);
vf=Properties(3);
% vg=Properties(4);
hfg=Properties(5);
cpf=Properties(6);
% cpg=Properties(7);
muf=Properties(8);
% mug=Properties(9);
kf=Properties(10);
% kg=Properties(11);
Prf=Properties(12);
% Prg=Properties(13);
% sigf=Properties(14);
% Betaf=Properties(15);
% Temp2=Properties(16);

end