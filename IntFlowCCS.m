function [Re,NuD,hbar] = IntFlowCCS(Fluid,mdot,D)
%INTFLOWCCS Summary of this function goes here
%   Detailed explanation goes here
 
Re= 4*mdot/(pi*D*Fluid.mu);

xfdh=10*D;
xfdt=xfdh; %#ok<*NASGU> 

f=(0.79*log(Re)-1.64)^-2;
NuD=(f/8*(Re-1000)*Fluid.Pr)/(1+12.7*sqrt(f/8)*(Fluid.Pr^(2/3)-1));
hbar=NuD*Fluid.k/D;
end

