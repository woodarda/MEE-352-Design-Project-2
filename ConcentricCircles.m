function [Run90CT,Run150CT] = ConcentricCircles(Run90CT,Run150CT,Chem,Water,Di,Do,D,ks)
%CONCENTRICCIRCLES Summary of this function goes here
%   Detailed explanation goes here
%% HEAT EXCHANGER: CONCENTRIC CIRCLES
% Calculates NTU if the flow is parallel
ParFlowNTU=@(Eff,Cr) -log(1-Eff*(1+Cr))/(1+Cr);

Run90CT.ParNTU=ParFlowNTU(Run90CT.Eff,Run90CT.Cr);
Run150CT.ParNTU=ParFlowNTU(Run150CT.Eff,Run150CT.Cr);

% Calculates NTU if the flow is Counter
CountFlowNTU=@(Eff,Cr) log((Eff-1)/(Eff*Cr-1))/(Cr-1);

% Runs CountNTU for both 90 and 150 scenarios
Run90CT.CountNTU=CountFlowNTU(Run90CT.Eff,Run90CT.Cr);
Run150CT.CountNTU=CountFlowNTU(Run150CT.Eff,Run150CT.Cr);

% Calculate Reynolds number, Nusselt number and convection coefficient 
% of the hot fluid
[Run90CT.Re_h,Run90CT.NuD_h,Run90CT.h_h,Run90CT.xfdh_h]=IntFlowCCS(Chem,Chem.mdot1,Di);
[Run150CT.Re_h,Run150CT.NuD_h,Run150CT.h_h,Run150CT.xfdh_h]=IntFlowCCS(Chem,Chem.mdot2,Di);

% Calculate Reynolds number, Nusselt number and convection coefficient 
% of the cold fluid
[Run90CT.Re_c,Run90CT.NuD_c,Run90CT.h_c,Run90CT.xfdh_c]=IntFlowCCS(Water,Water.mdot,D-Do);
[Run150CT.Re_c,Run150CT.NuD_c,Run150CT.h_c,Run150CT.xfdh_c]=IntFlowCCS(Water,Water.mdot,D-Do);

% Function to calculate Rconv
RConvFun=@(h,D) (h*pi*D)^-1;

% Calculate Req for the 90 min scenario
Run90CT.Rconvh=RConvFun(Run90CT.h_h,Di);
Run90CT.Rcond=log(Do/Di)/(2*pi*ks);
Run90CT.Rconvc=RConvFun(Run90CT.h_c,Do);

R=[Run90CT.Rconvh,Run90CT.Rcond,Run90CT.Rconvc];

Run90CT.Req=sum(R);

% Calculate Req for the 150 min scenario
Run150CT.Rconvh=RConvFun(Run150CT.h_h,Di);
Run150CT.Rcond=log(Do/Di)/(2*pi*ks);
Run150CT.Rconvc=RConvFun(Run150CT.h_c,Do);

R=[Run150CT.Rconvh,Run150CT.Rcond,Run150CT.Rconvc];

Run150CT.Req=sum(R);

% Func to calculate L
LFunc=@(Req,NTU,Cmin) Req*(NTU*Cmin);

% Calculate L for run 90 Counter Flow
Run90CT.LCount=LFunc(Run90CT.Req,Run90CT.CountNTU,Run90CT.Cmin);

% Calculate L for run 150 Parallel Flow
Run150CT.LPar=LFunc(Run150CT.Req,Run150CT.ParNTU,Run150CT.Cmin);
% Calculate L for run 150 Counter Flow
Run150CT.LCount=LFunc(Run150CT.Req,Run150CT.CountNTU,Run150CT.Cmin);
end

