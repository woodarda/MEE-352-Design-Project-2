function [Run90CF,Run150CF] = CrossFlow(Run90CF,Run150CF,Chem,Water,Di,Do,D,ks)
%CROSSFLOW Summary of this function goes here
%   Detailed explanation goes here

%% Calculate the NTU
Run90CF.Rho=Run90CF.v^-1;

[Run90CF.Cmin,Run90CF.Cmax,Run90CF.Cr,Run90CF.NTU]=NTU_CF(Run90CF.Eff,...
    Run90CF.Ch,Run90CF.Cc,Run90CF.CcMixed);

[Run150CF.Cmin,Run150CF.Cmax,Run150CF.Cr,Run150CF.NTU]=...
    NTU_CF(Run150CF.Eff, Run150CF.Ch,Run150CF.Cc,Run150CF.CcMixed);

% Calculate Reynolds number, Nusselt number and convection coefficient 
% of the hot fluid
[Run90CF.Re_h,Run90CF.NuD_h,Run90CF.h_h,Run90CF.xfdh_h]=IntFlowCCS(Chem,Chem.mdot1,Di);
[Run150CF.Re_h,Run150CF.NuD_h,Run150CF.h_h,Run150CF.xfdh_h]=IntFlowCCS(Chem,Chem.mdot2,Di);

% Calculate Reynolds number, Nusselt number and convection coefficient 
% of the cold fluid



end

