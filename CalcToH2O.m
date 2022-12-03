function [To_des] = CalcToH2O(To_des,Ti_des,mdot_des,mdot_known,cpf_known)
%CALCTOH2O Summary of this function goes here
%   Detailed explanation goes here
eps=100;
while eps>0.05
Tbar=(To_des+Ti_des)*0.5;
[Temp, P,vf, hfg,cpf_des,muf,kf,Prf]=AW_Interpolation(Tbar); %#ok<ASGLU> 

To_des=mdot_known*cpf_known/(mdot_des*cpf_des)+Ti_des;

Tbara=To_des/2+Ti_des/2;

eps=abs(Tbara-Tbar)/Tbara*100;
end

end

